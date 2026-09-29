"""The 1/q constant of κ for p = (q − 1)/2 (V30, V33).

q(κ((q−1)/2 / q) − κ₀) from the FFT route (`bulbford.taylor.kappa_fft`) for
q = 1025 … 8193, fitted by C + D/q + E/q². The phase alone contributes
`phase_curvature(q, HALF)` = i/π; the rest, R = C − i/π, is what the phase does
not explain. For p = 1 the phase contributes 1/(2πi), which is all of C (V28b).
V33 predicts R = iπ·dκ₀/dδ: the root sits at δ = iπ/q + … off c₀ = −3/4, and the
horn map of f + δ moves κ₀ by δ·dκ₀/dδ (`remainder`, at q = 10⁸ + 1). For p = 1,
δ = π²/q² + …, so the same term enters only at q⁻².
Writes experiments/data/kappa_constant.json.

    PYTHONPATH=kernel python3 experiments/scripts/kappa_constant.py
"""
from __future__ import annotations

import json

import mpmath as mp
import numpy as np

from bulbford.implosion import DPS, HALF, ONE, kappa0, kappa0_slope, phase_curvature, remainder
from bulbford.taylor import kappa_fft
from paths import DATA

OUT = DATA / "kappa_constant.json"
QS = (1025, 1537, 2049, 3073, 4097, 6145, 8193)


def fit(qs, ys, terms: int) -> tuple[complex, float]:
    X = np.array([[q ** -k for k in range(terms)] for q in qs])
    coef, *_ = np.linalg.lstsq(X, np.array(ys), rcond=None)
    return complex(coef[0]), float(np.max(np.abs(X @ coef - np.array(ys))))


def main() -> None:
    with mp.workdps(DPS):
        k0 = complex(kappa0(germ=HALF))
        slope = complex(kappa0_slope(HALF))
        big = 10**8 + 1
        predicted = complex(remainder((big - 1) // 2, big))
        curvature = {"p = 1": complex(phase_curvature(10**4, ONE)), "p = (q-1)/2": complex(phase_curvature(10**4, HALF))}
    rows = [{"q": q, "q_times_gap": q * (complex(kappa_fft((q - 1) // 2, q).coeffs[2]) - k0)} for q in QS]
    ys = [r["q_times_gap"] for r in rows]
    fits = {f"{t} terms": fit(QS, ys, t) for t in (2, 3)}
    C = fits["3 terms"][0]
    R = C - curvature["p = (q-1)/2"]
    for r in rows:
        print(f"q = {r['q']:5d}  q(κ − κ₀) = {r['q_times_gap']:.6f}")
    for k, (c, res) in fits.items():
        print(f"fit {k}: C = {c:.6f}  max residual {res:.1e}")
    print(f"phase curvature: p = 1 {curvature['p = 1']:.10f}, p = (q-1)/2 {curvature['p = (q-1)/2']:.10f}")
    print(f"R = C − i/π = {R:.6f}")
    print(f"dκ₀/dδ = {slope:.10f}, predicted R = iπ·dκ₀/dδ = {predicted:.10f}, gap {abs(R - predicted):.1e}")
    c = lambda z: [z.real, z.imag]  # noqa: E731
    OUT.write_text(json.dumps({
        "kappa0": c(k0), "curvature": {k: c(v) for k, v in curvature.items()},
        "rows": [{"q": r["q"], "q_times_gap": c(r["q_times_gap"])} for r in rows],
        "fits": {k: {"C": c(v[0]), "max_residual": v[1]} for k, v in fits.items()}, "R": c(R),
        "kappa0_slope": c(slope), "R_predicted": c(predicted),
    }, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
