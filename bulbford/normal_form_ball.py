"""P4 normal form of f(z) = λ₀z + z² in Arb ball arithmetic (python-flint): certified a, ι, κ in O(q² / B · M(q) + qB).

The recursion (see normal_form.py) needs the squares S_k = (H²)_k = Σ_{i+j=k} H_i H_j online: H_k depends on S_k, which
depends on H_1 … H_{k−1}.  Blocked (semi-relaxed) evaluation with block size B and frontier F (a multiple of B):

  * C[k] holds the ordered pairs with both indices < F.  When the block X = Σ_{F ≤ i < F+B} H_i z^i is complete,
    C += 2·P_F·X + X²  (P_F = Σ_{i<F} H_i z^i), one truncated series product in C (acb_poly_mullow).
  * For F ≤ k < F + B (and F ≥ B) the missing ordered pairs have one index in [F, k); the other is then < B, so
    S_k = C[k] + 2 Σ_{i=F}^{k−1} H_i H_{k−i}  — O(B) scalar operations.  The first block (F = 0) is summed directly.

All arithmetic is in balls, so the returned radii are rigorous error bounds."""
from __future__ import annotations
from typing import NamedTuple
from flint import acb, acb_series, arb, ctx


class BallNormalForm(NamedTuple):
    """Resonant coefficients and derived invariants, all evaluated at the working precision (balls)."""
    p: int
    q: int
    b1: acb
    b2: acb
    a: acb          # [z^{q+1}] f^q = q·b₁
    iota: acb       # ι_{p/q} = (q²−1)/(2q) + b₂/(q b₁²)
    kappa: acb      # κ(p/q) = (ι − ½)/q


def _root_power(p: int, q: int, k: int) -> acb:
    """λ₀^k = exp(2πi·(pk mod q)/q), evaluated directly (no radius growth from repeated products)."""
    return (acb(0, 2) * arb.pi() * ((p * k) % q) / q).exp()


def _normal_form_ball(p: int, q: int, prec: int, B: int) -> BallNormalForm:
    L = 2 * q + 1
    old_prec, old_cap = ctx.prec, ctx.cap
    try:
        ctx.prec, ctx.cap = prec, L + 1
        lam = _root_power(p, q, 1)
        H = [acb(0)] * (L + 1)
        H[1] = acb(1)
        C = [acb(0)] * (L + 1)
        F = 0
        beta1 = beta2 = acb(0)
        for k in range(2, L + 1):
            if k >= F + B:                                   # fold the completed block [F, F+B) into C
                n = L - F + 1                                # C[F + i] is needed for i < n only
                X = acb_series(H[F:F + B], prec=n)           # X/z^F
                if F:
                    PX = (acb_series(H[:F], prec=n) * X).coeffs()        # (P_F · X)/z^F, truncated
                    for i, c in enumerate(PX):
                        C[F + i] += 2 * c
                if 2 * F <= L:
                    XX = (X * X).coeffs()                                # X²/z^{2F}
                    for i, c in enumerate(XX[:L - 2 * F + 1]):
                        C[2 * F + i] += c
                F += B
            if F == 0:
                S = sum((H[i] * H[k - i] for i in range(1, k)), acb(0))
            else:
                S = C[k] + 2 * sum((H[i] * H[k - i] for i in range(F, k)), acb(0))
            if k == q + 1:
                beta1 = S; continue                          # H_{q+1} = 0
            if k == 2 * q + 1:
                beta2 = S; continue                          # H_{2q+1} = 0
            rhs = -S
            if k > q + 1:
                m = k - q
                rhs += m * _root_power(p, q, m - 1) * beta1 * H[m]
            H[k] = rhs / (lam - _root_power(p, q, k))
        b1, b2 = beta1 / lam, beta2 / lam
        iota = acb(q * q - 1) / (2 * q) + b2 / (q * b1 ** 2)
        return BallNormalForm(p, q, b1, b2, q * b1, iota, (iota - acb(1) / 2) / q)
    finally:
        ctx.prec, ctx.cap = old_prec, old_cap


def normal_form_ball(p: int, q: int, prec: int = 256, B: int | None = None,
                     target_rad: float | None = 1e-40, max_prec: int = 1 << 15) -> BallNormalForm:
    """Certified normal form.  Radii are rigorous but pessimistic (worst-case accumulation through the long
    convolutions); if the radius of κ exceeds `target_rad`, the working precision is doubled and the computation
    repeated (up to `max_prec` bits).  `target_rad=None` returns the first result."""
    B = B or max(32, int(2 * (2 * q + 1) ** 0.5))
    while True:
        nf = _normal_form_ball(p, q, prec, B)
        if target_rad is None or float(nf.kappa.rad()) <= target_rad or prec >= max_prec:
            return nf
        prec *= 2
