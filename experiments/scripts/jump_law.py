"""Next move 2: the jump law of Im κ and the singularity profile of Re κ at several q (V33).

For each q (kappa_q<q>.json must exist):
  1. Jump law.  π m S_m = Σ_{q'} J_{q'} c_{q'}(m) + b/m² on windows [lo, hi]; J·q'² for q' ≤ 10 and the
     log-log slope of |J_{q'}| on q' = 3 … 10 (−2 is the Ford weight; drift shows as a different slope).
  2. Profile of Re κ.  C_m = g(m) Σ D_{q'} c_{q'}(m) + b/m⁴, weighted by m, for the profiles g of
     spectral.py. (A free Hölder-exponent scan is not used: at high m the flat sampling-noise floor
     favours the flattest profile.) The same m-windows at every q: a q-independent log amplitude D means
     a true log singularity, D rising with q in the high windows means one smoothed at the scale 1/q.
  3. Hilbert pair.  If κ were the boundary value of a function holomorphic above the real line, each
     Im-jump J would carry the Re-singularity −(J/π) log|δ|, i.e. the log-profile amplitude D = J/π.
  4. Real space.  Re κ against n smooth cosine modes plus one singular regressor: the truncated Brjuno sum
     B_T(x*) (continued-fraction organised), or a log sum Σ_{q' ≤ 20} w(q') Σ_{p'} log|2 sin π(x̃ − p'/q')| with
     fixed weights w = q'^{-1}, q'^{-3/2}, q'^{-2} (every rational, no continued-fraction organisation).
Writes experiments/data/jump_law.json.

    PYTHONPATH=kernel python3 experiments/scripts/jump_law.py 1009 2003
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import gcd

import numpy as np

from paths import DATA
from spectral import brjuno, coefficients, kappa_table, ramanujan_fit

QMAX = 10
PROFILES = {"log |δ| (1/m)": lambda m: 1.0 / m, "V-cusp (1/m²)": lambda m: 1.0 / m**2,
            "|δ| log|δ| (log m/m²)": lambda m: np.log(m) / m**2}
WINDOWS = ((4, 60), (8, 120), (12, 200), (20, 300))


def weighted_fit(y, lo, hi, g, weight=lambda m: m * 1.0):
    """ramanujan_fit on w·y with basis w·g·c_{q'}, background w/m⁴ (same weights for every profile)."""
    m = np.arange(1, len(y) + 1)
    return ramanujan_fit(weight(m) * y, lo, hi, QMAX, g=lambda mm: weight(mm) * g(mm),
                         background=lambda mm: weight(mm) / mm**4)


LOG_WEIGHTS = {"q'^-1": lambda b: 1.0 / b, "q'^-1.5": lambda b: b ** -1.5, "q'^-2": lambda b: 1.0 / b**2}


def log_sum(x: np.ndarray, w, qmax: int = 20) -> np.ndarray:
    return sum(w(b) * np.log(np.abs(2 * np.sin(np.pi * (x - a / b)))) for b in range(1, qmax + 1)
               for a in range(b) if gcd(a, b) == 1)


def real_space(q: int, bounded: int = 4) -> list[dict]:
    """Residual sd of Re κ after n cosine modes, and after adding each singular regressor."""
    p, xt, k = kappa_table(q)
    keep = np.minimum(p, q - p) > bounded
    x, y = xt[keep], k.real[keep]
    xs = [Fraction(min(a, q - a), q) for a in np.rint(x * q).astype(int)]     # x* = |x̃|: Re κ is even in x̃
    regressors = {**{f"Brjuno T={T}": np.array([brjuno(v, T) for v in xs]) for T in (5, 45, 300)},
                  **{f"log sum {name}": log_sum(x, w) for name, w in LOG_WEIGHTS.items()}}
    rows = []
    for n in (8, 16, 32):
        base = np.column_stack([np.ones_like(x)] + [np.cos(2 * np.pi * j * x) for j in range(1, n)])
        sd0 = float(np.std(y - base @ np.linalg.lstsq(base, y, rcond=None)[0]))
        row = {"modes": n, "sd": sd0}
        for name, r in regressors.items():
            A = np.column_stack([base, r]); c = np.linalg.lstsq(A, y, rcond=None)[0]
            row[name] = {"sd": float(np.std(y - A @ c)), "coef": float(c[-1])}
        rows.append(row)
    return rows


def analyse(q: int) -> dict:
    co = coefficients(q, q // 2)
    m, S, C = co["m"], co["S"], co["C"]
    tail = slice(q // 4, q // 2)
    out = {"q": q, "noise_rms": {"S": float(np.sqrt(np.mean(S[tail] ** 2))), "C": float(np.sqrt(np.mean(C[tail] ** 2)))},
           "jump": [], "profiles": [], "hilbert": []}
    for lo, hi in WINDOWS:
        J, r2 = ramanujan_fit(np.pi * m * S, lo, hi, QMAX)
        slope = float(np.polyfit(np.log(np.arange(3, QMAX + 1)), np.log(np.abs(J[2:])), 1)[0])
        out["jump"].append({"window": [lo, hi], "r2": r2, "J_q2": [float(J[b - 1] * b * b) for b in range(1, QMAX + 1)],
                            "loglog_slope_q3_10": slope})
        out["profiles"].append({"window": [lo, hi], **{name: weighted_fit(C, lo, hi, g)[1] for name, g in PROFILES.items()}})
        D, _ = weighted_fit(C, lo, hi, PROFILES["log |δ| (1/m)"])
        out["hilbert"].append({"window": [lo, hi], "D_q": [float(D[b - 1]) for b in range(1, QMAX + 1)],
                               "D_over_J_by_pi": [float(D[b - 1] / (J[b - 1] / np.pi)) for b in range(2, 8)],
                               "signal_at_hi": {"S": float(abs(S[hi - 1])), "C": float(abs(C[hi - 1]))}})
    out["real_space"] = real_space(q)
    return out


def report(r: dict) -> None:
    print(f"\n=== q = {r['q']} ===  noise rms (m ∈ [q/4, q/2]): S {r['noise_rms']['S']:.1e}, C {r['noise_rms']['C']:.1e}")
    for j in r["jump"]:
        print(f"J·q'² {j['window']} R²={j['r2']:.3f} slope(q'=3..10)={j['loglog_slope_q3_10']:+.2f}:",
              " ".join(f"{v:+.3f}" for v in j["J_q2"]))
    for pr, h in zip(r["profiles"], r["hilbert"]):
        print(f"Re κ {pr['window']}:", "  ".join(f"{k}: R²={v:.3f}" for k, v in pr.items() if k != "window"),
              "| D·q', q'=2..7:", " ".join(f"{h['D_q'][b-1]*b:+.4f}" for b in range(2, 8)),
              "| D/(J/π):", " ".join(f"{v:+.2f}" for v in h["D_over_J_by_pi"]),
              f"| |C_hi| {h['signal_at_hi']['C']:.1e}")
    for row in r["real_space"]:
        print(f"Re κ real space, {row['modes']:2d} modes: sd {row['sd']:.2e} |",
              " | ".join(f"{k}: {1 - v['sd'] / row['sd']:+.0%} sd ({1 - (v['sd'] / row['sd'])**2:+.0%} var), coef {v['coef']:+.4f}"
                         for k, v in row.items() if isinstance(v, dict)))


if __name__ == "__main__":
    results = [analyse(int(q)) for q in sys.argv[1:]]
    for r in results:
        report(r)
    (DATA / "jump_law.json").write_text(json.dumps(results, indent=1))
