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
