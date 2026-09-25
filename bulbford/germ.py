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
