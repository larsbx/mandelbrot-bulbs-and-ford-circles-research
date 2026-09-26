"""Least-squares power fits y(q) = Σ_{j≤J} c_j q^{−j} in mpmath (used for all q → ∞ extrapolations).
QR least squares (not normal equations): the columns q^{−j} are badly scaled, cond(AᵀA) ~ q_max^{2J}."""
from __future__ import annotations
import mpmath as mp


def power_fit(rows, J: int, dps: int = 80):
    """rows: iterable of (q, y) with y real or complex.  Returns [c_0, …, c_J] (mpmath numbers)."""
    rows = sorted(rows)
    with mp.workdps(dps):
        A = mp.matrix([[mp.mpc(mp.mpf(q) ** (-j)) for j in range(J + 1)] for q, _ in rows])
        y = mp.matrix([mp.mpc(mp.mpmathify(v)) for _, v in rows])
        x, _ = mp.qr_solve(A, y)
        return [x[j] for j in range(J + 1)]


def linear_limit(xs):
    """Limit of a sequence x_k ∈ ℂ whose differences obey Δ_{k+1} = M Δ_k with M real-linear on ℂ ≅ ℝ² — this covers
    holomorphic (Δ ↦ νΔ) and antilinear (Δ ↦ νΔ̄, eigenvalues ±|ν|) contractions alike, where per-component Aitken
    is not exact.  M is fitted from the last three differences.  Returns (limit, eigenvalues of M)."""
    vec = lambda z: mp.matrix([mp.re(z), mp.im(z)])
    d = [vec(b - a) for a, b in zip(xs[-4:], xs[-3:])]
    A, B = mp.matrix(2, 2), mp.matrix(2, 2)
    for j in range(2):
        A[:, j], B[:, j] = d[j], d[j + 1]
    M = B * A ** -1
    rem = (mp.eye(2) - M) ** -1 * M * d[2]
    return xs[-1] + mp.mpc(rem[0], rem[1]), mp.eig(M)[0]
