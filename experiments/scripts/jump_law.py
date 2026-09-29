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
Writes experiments/data/jump_law.json.

    PYTHONPATH=kernel python3 experiments/scripts/jump_law.py 1009 2003
"""
from __future__ import annotations

import json
import sys

import numpy as np

from paths import DATA
from spectral import coefficients, ramanujan_fit

QMAX = 10
PROFILES = {"log |δ| (1/m)": lambda m: 1.0 / m, "V-cusp (1/m²)": lambda m: 1.0 / m**2,
            "|δ| log|δ| (log m/m²)": lambda m: np.log(m) / m**2}
WINDOWS = ((4, 60), (8, 120), (12, 200), (20, 300))


def weighted_fit(y, lo, hi, g, weight=lambda m: m * 1.0):
    """ramanujan_fit on w·y with basis w·g·c_{q'}, background w/m⁴ (same weights for every profile)."""
    m = np.arange(1, len(y) + 1)
    return ramanujan_fit(weight(m) * y, lo, hi, QMAX, g=lambda mm: weight(mm) * g(mm),
                         background=lambda mm: weight(mm) / mm**4)


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


if __name__ == "__main__":
    results = [analyse(int(q)) for q in sys.argv[1:]]
    for r in results:
        report(r)
    (DATA / "jump_law.json").write_text(json.dumps(results, indent=1))
