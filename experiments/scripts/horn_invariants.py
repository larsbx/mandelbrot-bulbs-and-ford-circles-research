"""The horn maps against the two V25 families (V28 for p = 1, V29 for p = (q − 1)/2).

For each family's germ (`bulbford.implosion.ONE`, `HALF`) this finds, for each
named cycle of V27, the fixed point of the Lavaurs map L_{σ∞} seeded by the cycle.
It tabulates |1 − L'| against |1 − F'_q| and the gap |F'_q − L'| for three values
of q. It computes κ₀ = end·a₂/(2πi a₁²) from the horn map's Fourier coefficients
on two lines, and tabulates q(κ(p/q) − κ₀) from the exact index and
experiments/data/kappa_q1009.json. Writes experiments/data/horn_invariants.json.

    PYTHONPATH=kernel python3 experiments/scripts/horn_invariants.py
"""
from __future__ import annotations

import json

import mpmath as mp
import numpy as np

from bulbford.cycles import cycle_points
from bulbford.implosion import DPS, HALF, ONE, fixed_point_of_cycle, kappa0, phase, phi_in
from bulbford.index import kappa
from multiplier_growth import named_angles
from paths import DATA

OUT = DATA / "horn_invariants.json"

FAMILIES = (
    {"germ": ONE, "p": lambda q: 1, "cycles_q": (64, 128, 256), "kappa_q": (64, 128, 256, 512), "lines": (2.8, 3.2),
     "seeds": 4},
    {"germ": HALF, "p": lambda q: (q - 1) // 2, "cycles_q": (513, 1025, 2049), "kappa_q": (65, 129, 257, 513),
     "lines": (-2.4, -2.9), "seeds": 3},
)


def c(x) -> list[float]:
    x = complex(x)
    return [x.real, x.imag]


def cycles(family) -> list[dict]:
    germ, rows = family["germ"], []
    for q in family["cycles_q"]:
        p = family["p"](q)
        for name, angle in named_angles(p, q).items():
            z = cycle_points(p, q, angle)[1]
            F = complex(np.prod(2 * z))
            found = fixed_point_of_cycle(z, germ=germ, seeds=family["seeds"])
            w, slope = found[0]
            rows.append({"q": q, "cycle": name, "seeds_agreeing": sum(abs(s - slope) < 1e-20 for _, s in found),
                         "seeds_converged": len(found), "Z": c(phi_in(w, germ)[0] + phase(germ)),
                         "one_minus_L": float(abs(1 - slope)), "one_minus_F": abs(1 - F), "gap": float(abs(slope - F))})
            print(f"{germ.name}: q = {q:4d} {name:13s} |1 − L'| = {rows[-1]['one_minus_L']:.9f}  |1 − F'| = {abs(1 - F):.6f}"
                  f"  q·gap = {q * rows[-1]['gap']:.1f}  q²·gap = {q * q * rows[-1]['gap']:.0f}", flush=True)
    return rows


def kappas(family, k0) -> list[dict]:
    p_of = family["p"]
    rows = [{"q": q, "p": p_of(q), "kappa": c(kappa(p_of(q), q))} for q in family["kappa_q"]]
    table = json.loads((DATA / "kappa_q1009.json").read_text())
    rows.append({"q": 1009, "p": p_of(1009), "kappa": next(r["kappa"] for r in table if r["p"] == p_of(1009))})
    for r in rows:
        r["q_times_gap"] = c(r["q"] * (complex(*r["kappa"]) - complex(k0)))
    for a, b in zip(rows, rows[1:]):
        b["richardson"] = c((b["q"] * complex(*b["q_times_gap"]) - a["q"] * complex(*a["q_times_gap"])) / (b["q"] - a["q"]))
    for r in rows:
        print(f"{family['germ'].name}: q = {r['q']:4d}  q(κ − κ₀) = {complex(*r['q_times_gap']):.6f}"
              + (f"  Richardson {complex(*r['richardson']):.6f}" if "richardson" in r else ""), flush=True)
    return rows


def main() -> None:
    out = []
    with mp.workdps(DPS):
        for family in FAMILIES:
            germ = family["germ"]
            k0, k1 = (kappa0(y, germ) for y in family["lines"])
            print(f"{germ.name}: σ∞ = {mp.nstr(phase(germ), 12)}, κ₀ = {mp.nstr(k0, 15)} (lines differ by {mp.nstr(abs(k0 - k1), 3)})")
            out.append({"family": germ.name, "phase": c(phase(germ)), "kappa0": c(k0), "kappa0_line_drift": float(abs(k0 - k1)),
                        "cycles": cycles(family), "kappa": kappas(family, k0)})
    OUT.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
