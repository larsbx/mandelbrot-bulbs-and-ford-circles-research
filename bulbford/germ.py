"""Truncated power series (mpmath lists, index = degree) and holomorphic indices of iterates of germs."""
from __future__ import annotations
import mpmath as mp


def mul(a, b, L):
    return [mp.fsum(a[i] * b[k - i] for i in range(max(0, k - len(b) + 1), min(k, len(a) - 1) + 1)) for k in range(L + 1)]


def compose(f, g, L):
    """f(g(W)) for g(0) = 0, degrees ≤ L (Horner)."""
    out = [mp.mpc(0)] * (L + 1)
    for c in reversed(f[:L + 1]):
        out = mul(out, g, L); out[0] += c
    return out


def exp_series(a, L):
    """exp(Σ a_k W^k), a_0 = 0, via E' = A'E."""
    e = [mp.mpc(1)] + [mp.mpc(0)] * L
    for k in range(1, L + 1):
        e[k] = mp.fsum(j * a[j] * e[k - j] for j in range(1, min(k, len(a) - 1) + 1)) / k
    return e


def horn_germ(a, L, mu=1):
    """G(W) = μ·W·exp(2πi Σ_{n≥1} a_n W^n) from upper horn-map coefficients a (a₀ must vanish)."""
    ex = exp_series([mp.mpc(0)] + [2j * mp.pi * x for x in a[1:]], L - 1)
    return [mp.mpc(0)] + [mu * c for c in ex[:L]]


def iterate_index(G, k):
    """ι of G^k at 0 when G'(0) is a primitive k-th root of unity: W − G^k(W) = W^{k+1}P(W), ι = [W^k] 1/P."""
    L = 2 * k + 1
    Gk = [mp.mpc(0), mp.mpc(1)] + [mp.mpc(0)] * (L - 1)
    for _ in range(k):
        Gk = compose(G, Gk, L)
    D = [(1 if d == 1 else 0) - Gk[d] for d in range(L + 1)]
    P = D[k + 1:]
    inv = [1 / P[0]] + [mp.mpc(0)] * k
    for m in range(1, k + 1):
        inv[m] = -mp.fsum(P[j] * inv[m - j] for j in range(1, m + 1)) / P[0]
    return inv[k], max((abs(x) for x in D[2:k + 1]), default=mp.mpf(0))


def normal_form_index(g, p):
    """ι(G^p) for a germ G(z) = Σ_{j≥1} g_j z^j (g[0] = 0) whose multiplier μ = g[1] is a primitive p-th root of
    unity, via the equivariant normal form G∘H = H∘N, N = μz(1 + b₁z^p + b₂z^{2p}) (P4 for general germs):
      H_k(μ − μ^k) = [m = k−p ≥ 1]·m μ^m b₁ H_m + [m = k−2p ≥ 1]·μ^m(m b₂ + C(m,2) b₁²) H_m − Σ_{j≥2} g_j [H^j]_k,
    b₁ from k = p+1, b₂ from k = 2p+1 (H_{p+1} = H_{2p+1} = 0);  ι(N^p) = (p²−1)/(2p) + b₂/(p b₁²).
    Cost O(p³) scalar operations (powers of H kept in a table), versus O(p⁴) for composing series p times."""
    L = 2 * p + 1
    mu = g[1]
    gj = lambda j: g[j] if j < len(g) else 0
    H = [mp.mpc(0)] * (L + 1); H[1] = mp.mpc(1)
    P = [[mp.mpc(0)] * (L + 1) for _ in range(L + 1)]     # P[j][m] = [H^j]_m
    P[1][1] = mp.mpc(1)
    b1 = b2 = mp.mpc(0)
    for k in range(2, L + 1):
        for j in range(2, k + 1):
            P[j][k] = mp.fsum(H[i] * P[j - 1][k - i] for i in range(1, k - j + 2))
        S = mp.fsum(gj(j) * P[j][k] for j in range(2, k + 1))
        if k == p + 1:
            b1 = S / mu; continue                          # H_{p+1} = 0
        if k == 2 * p + 1:
            b2 = S / mu; continue                          # H_{2p+1} = 0 (m = p+1 term vanishes)
        rhs = -S
        m = k - p
        if m >= 1:
            rhs += m * mu ** m * b1 * H[m]
        m = k - 2 * p
        if m >= 1:
            rhs += mu ** m * (m * b2 + mp.binomial(m, 2) * b1 ** 2) * H[m]
        H[k] = rhs / (mu - mu ** k)
        P[1][k] = H[k]
    return mp.mpf(p * p - 1) / (2 * p) + b2 / (p * b1 ** 2)
