"""The horn map of w + w² against the p = 1 roots (V28).

For the named cycles of V27 (p = 1) this finds the fixed point of the Lavaurs map
L_{iπ} seeded by the cycle, and tabulates |1 − L'| against |1 − F'_q| and the gap
|F'_q − L'| for q = 64, 128, 256. It computes κ₀ = a₂/(2πi a₁²) from the Fourier
coefficients of the horn map on two lines, and tabulates q(κ(1/q) − κ₀) from the
exact index (q ≤ 512) and experiments/data/kappa_q1009.json. Writes experiments/data/horn_invariants.json.

    PYTHONPATH=kernel python3 experiments/scripts/horn_invariants.py
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import numpy as np

from bulbford.cycles import cycle_points
from bulbford.implosion import DPS, fixed_point_of_cycle, kappa0, phi_in
from bulbford.index import kappa
from multiplier_growth import named_angles
from paths import DATA

OUT = DATA / "horn_invariants.json"


def c(x) -> list[float]:
    x = complex(x)
    return [x.real, x.imag]


def cycles() -> list[dict]:
    rows = []
    for q in (64, 128, 256):
        for name, angle in named_angles(1, q).items():
            z = cycle_points(1, q, angle)[1]
            F = complex(np.prod(2 * z))
            found = fixed_point_of_cycle(z)
            w, slope = found[0]
            rows.append({"q": q, "cycle": name, "seeds_agreeing": sum(abs(s - slope) < 1e-20 for _, s in found),
                         "Z": c(phi_in(w)[0] + 1j * mp.pi), "one_minus_L": float(abs(1 - slope)),
                         "one_minus_F": abs(1 - F), "gap": float(abs(slope - F))})
            print(f"q = {q:3d} {name:13s} |1 − L'| = {rows[-1]['one_minus_L']:.9f}  |1 − F'| = {abs(1 - F):.6f}"
                  f"  q²·gap = {q * q * rows[-1]['gap']:.0f}", flush=True)
    return rows


def kappas(k0) -> list[dict]:
    rows = [{"q": q, "kappa": c(kappa(1, q))} for q in (64, 128, 256, 512)]
    table = json.loads((DATA / "kappa_q1009.json").read_text())
    rows.append({"q": 1009, "kappa": next(r["kappa"] for r in table if r["p"] == 1)})
    for r in rows:
        r["q_times_gap"] = c(r["q"] * (complex(*r["kappa"]) - complex(k0)))
        print(f"q = {r['q']:4d}  q(κ − κ₀) = {complex(*r['q_times_gap']):.6f}", flush=True)
    return rows


def main() -> None:
    with mp.workdps(DPS):
        k0 = kappa0(2.8)
        drift = abs(k0 - kappa0(3.2))
        print(f"κ₀ = {mp.nstr(k0, 15)}  (lines Im Z = 2.8, 3.2 differ by {mp.nstr(drift, 3)})")
        out = {"kappa0": c(k0), "kappa0_line_drift": float(drift), "one_over_2pi": float(1 / (2 * mp.pi)),
               "cycles": cycles(), "kappa": kappas(k0)}
    OUT.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
