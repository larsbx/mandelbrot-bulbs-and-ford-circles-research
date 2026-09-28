"""Exploratory spike: bridges from the bulb-ford study to arithmetic (docs/bridges-spike.md).

B1  The Fourier spectrum of κ at q = 1009 in the variable p̄/q against the variable p/q.
B2  Ramanujan-sum structure of that spectrum: π m S_m = Σ_{q'} J_{q'} c_{q'}(m) recovers the Im κ
    jump law at every denominator at once; the same model for Re κ with a slope-jump (V-cusp) basis.
B3  Norms of the parabolic coefficient a_q = [z^{q+1}] f^q in ℤ[ζ_q]: primes ℓ > 2q + 2 with
    ord_ℓ(λ) = q for λ ∈ {2, 4, −2} divide N(a_q), once per such λ (the maps with those
    multipliers are w², w² − 2 and w² − 2, which are exactly solvable).
B4  Dedekind sums s(p, q) against the G residual left by a function of x*.

Reads experiments/data/kappa_q1009.json; writes experiments/data/bridges_spike.json.

    PYTHONPATH=kernel python3 experiments/scripts/bridges_spike.py
"""
from __future__ import annotations

import json
from fractions import Fraction
from math import gcd

import numpy as np
from sympy import cyclotomic_poly, factorint, n_order

from bulbford.norms import norm_a
from paths import DATA

Q = 1009
OUT = DATA / "bridges_spike.json"
V21 = {2: -0.034, 3: -0.0186, 4: -0.0105, 5: -0.0049, 6: -0.0029}


def kappa_table() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(p, p̄/q, κ) over all units p mod Q, from the p ≤ Q/2 table and κ(Q − p) = conj κ(p)."""
    rows = json.loads((DATA / "kappa_q1009.json").read_text())
    half = {r["p"]: complex(*r["kappa"]) for r in rows}
    full = {**half, **{Q - p: k.conjugate() for p, k in half.items()}}
    p = np.array(sorted(full))
    return p, np.array([pow(int(a), -1, Q) for a in p]) / Q, np.array([full[a] for a in p])


def fourier(x: np.ndarray, y: np.ndarray, m: np.ndarray, trig) -> np.ndarray:
    return np.array([2 * np.mean(y * trig(2 * np.pi * k * x)) for k in m])


def b1(p, xt, kappa) -> dict:
    m = np.arange(1, 7)
    centred = kappa - kappa.mean()
    amp = lambda x: np.abs([np.mean(centred * np.exp(-2j * np.pi * k * x)) for k in m])
    in_pbar, in_p = amp(xt), amp(p / Q)
    return {"modes": m.tolist(), "abs_c_in_pbar": in_pbar.tolist(), "abs_c_in_p": in_p.tolist(),
            "ratio_l2": float(np.linalg.norm(in_pbar) / np.linalg.norm(in_p))}


def ramanujan(qq: int, m: np.ndarray) -> np.ndarray:
    return sum(np.cos(2 * np.pi * a * m / qq) for a in range(qq) if gcd(a, qq) == 1)


def ramanujan_fit(T: np.ndarray, lo: int, hi: int, qmax: int) -> tuple[np.ndarray, float]:
    """Least squares T(m) = Σ_{q' ≤ qmax} X_{q'} c_{q'}(m) + b/m² on m ∈ [lo, hi]."""
    m = np.arange(lo, hi + 1)
    X = np.column_stack([ramanujan(b, m) for b in range(1, qmax + 1)] + [1.0 / m**2])
    c, *_ = np.linalg.lstsq(X, T[lo - 1:hi], rcond=None)
    r = T[lo - 1:hi] - X @ c
    return c[:qmax], float(1 - r.var() / T[lo - 1:hi].var())


WINDOWS = ((4, 60, 10), (8, 120, 12), (12, 200, 14), (20, 300, 16))


def b2(p, xt, kappa) -> dict:
    keep = np.minimum(p, Q - p) > 4                       # the bounded-p term is separate (C1‴)
    x, k = xt[keep], kappa[keep]
    m = np.arange(1, WINDOWS[-1][1] + 1)
    S = fourier(x, k.imag, m, np.sin)
    C = fourier(x, k.real - k.real.mean(), m, np.cos)
    fits = lambda T, power: [
        {"window": [lo, hi], "qmax": qm, "r2": r2, "scaled": [float(c[b - 1] * b**power) for b in range(1, 9)]}
        for lo, hi, qm in WINDOWS for c, r2 in [ramanujan_fit(T, lo, hi, qm)]]
    return {"jump_J_times_q2": fits(np.pi * m * S, 2),
            "cusp_D_times_q": fits(2 * np.pi**2 * m**2 * C, 1),
            "V21_J_times_q2": {b: V21[b] * b * b for b in V21},
            "S_times_m_first12": (S[:12] * m[:12]).tolist()}


def b3(qs=range(3, 21)) -> dict:
    rows = []
    for q in qs:
        n = norm_a(q)
        primes = {l for base in (2, 4, -2) for l in factorint(abs(int(cyclotomic_poly(q, base)))) if l > 2 * q + 2}
        for l in sorted(primes):
            lams = sorted({x % l for x in (2, 4, -2) if n_order(x % l, l) == q})
            v = 0
            while n % l**(v + 1) == 0:
                v += 1
            rows.append({"q": q, "l": l, "lambdas": lams, "predicted": len(lams), "valuation": v})
    return {"rows": rows, "all_hold": all(r["valuation"] >= r["predicted"] for r in rows),
            "all_equal": all(r["valuation"] == r["predicted"] for r in rows)}


def dedekind(h: int, k: int) -> Fraction:
    """s(h, k) by the reciprocity algorithm."""
    s, sign, h = Fraction(0), 1, h % k
    while k > 1 and h:
        s += sign * (Fraction(h * h + k * k + 1, 12 * h * k) - Fraction(1, 4))
        h, k, sign = k % h, h, -sign
    return s


def b4() -> dict:
    rows = json.loads((DATA / "kappa_q1009.json").read_text())
    p = np.array([r["p"] for r in rows]); G = np.array([r["G"] for r in rows])
    xs = np.array([min(pow(int(a), -1, Q), Q - pow(int(a), -1, Q)) for a in p]) / Q
    bins = np.minimum((xs * 200).astype(int), 99)
    res = G - np.array([G[bins == b].mean() for b in bins])
    s = np.array([float(dedekind(int(a), Q)) for a in p])
    sd_after = lambda *cols: float(np.std(res - np.column_stack(cols + (np.ones_like(G),)) @
                                          np.linalg.lstsq(np.column_stack(cols + (np.ones_like(G),)), res, rcond=None)[0]))
    return {"residual_sd": float(res.std()), "corr_s": float(np.corrcoef(res, s)[0, 1]),
            "sd_after_inv_p": sd_after(1.0 / p), "sd_after_s": sd_after(s / Q), "sd_after_both": sd_after(1.0 / p, s / Q)}


def main() -> None:
    p, xt, kappa = kappa_table()
    out = {"B1": b1(p, xt, kappa), "B2": b2(p, xt, kappa), "B3": b3(), "B4": b4()}
    OUT.write_text(json.dumps(out, indent=1))
    print(f"B1 ‖c‖ in p̄/q over p/q (m ≤ 6): {out['B1']['ratio_l2']:.1f}")
    for f in out["B2"]["jump_J_times_q2"]:
        print(f"B2 J·q'² {f['window']} R²={f['r2']:.3f}:", " ".join(f"{v:+.3f}" for v in f["scaled"]))
    print("B2 cusp-model R²:", [round(f["r2"], 3) for f in out["B2"]["cusp_D_times_q"]])
    print(f"B3 {len(out['B3']['rows'])} (q, ℓ) checks, q ≤ 20: hold={out['B3']['all_hold']} equality={out['B3']['all_equal']}")
    print("B4", {k: round(v, 5) for k, v in out["B4"].items()})


if __name__ == "__main__":
    main()
