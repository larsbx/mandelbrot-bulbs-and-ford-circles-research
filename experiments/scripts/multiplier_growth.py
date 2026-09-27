"""How the terms of the named cycles grow with q (V27).

A cycle of period k contributes −k/(1 − F') with F' = ρ^{q/k}. For the two V25
families (p = 1 and p = (q − 1)/2) this tabulates |1 − F'| for the flank cycle,
the two ring-2 cycles and the widest-arc intruders, from single verified orbits
(`bulbford.cycles.cycle_through`), and fits the model a + b/q + c/q² on the
nodes q ≥ 64 (a fitted statistic, not a limit). The β fixed point (angle 0) has
the exact multiplier (2 − λ₀)^q and is tabulated from that closed form.
Writes experiments/data/multiplier_growth.json.

    PYTHONPATH=kernel python3 experiments/scripts/multiplier_growth.py [Q]
"""
from __future__ import annotations

import cmath
import json
import math
import sys
from pathlib import Path

import numpy as np

from bulbford.cycles import cycle_through
from bulbford.wake import mechanical
from paths import DATA

OUT = DATA / "multiplier_growth.json"


def named_angles(p: int, q: int) -> dict[str, int]:
    M, w = 2**q - 1, p - 1
    A = [mechanical(p, q, r) for r in range(q)]
    named = {"flank": (A[(w - 1) % q] + 1) % M, "ring2_lower": (A[(w - 2) % q] + 1) % M,
             "ring2_upper": (A[(w + 3) % q] - 1) % M}
    if p == 1:
        named["intruder_U-2"] = (A[(w - 1) % q] - 1) % M  # upper end of arc j = −2 (width 2^{q−2}/M)
    else:
        named["intruder_L1"] = (A[w + 1] + 1) % M  # lower end of arc j = +1 (width 2^{q−2}/M when p̄ = q − 2)
    return named


def beta_distance(p: int, q: int) -> float:
    return abs(1 - (2 - cmath.exp(2j * math.pi * p / q)) ** q)


def fit(qs: list[int], ys: list[float]) -> dict:
    X = np.array([[1.0, 1.0 / q, 1.0 / q**2] for q in qs])
    coef, *_ = np.linalg.lstsq(X, np.array(ys), rcond=None)
    residual = float(np.max(np.abs(X @ coef - np.array(ys))))
    return {"a": float(coef[0]), "b": float(coef[1]), "c": float(coef[2]), "nodes": qs, "max_residual": residual}


def family(name: str, p_of_q, qs: list[int]) -> dict:
    rows = []
    for q in qs:
        p = p_of_q(q)
        row = {"q": q, "p": p, "beta": beta_distance(p, q)}
        for label, angle in named_angles(p, q).items():
            t = cycle_through(p, q, angle)
            row[label] = abs(1 - t.multiplier ** (q // t.period))
        rows.append(row)
    labels = [k for k in rows[0] if k not in ("q", "p", "beta")]
    nodes = [r for r in rows if r["q"] >= 64]
    fits = {k: fit([r["q"] for r in nodes], [r[k] for r in nodes]) for k in labels}
    if all(r["p"] == 1 for r in rows):  # for p = (q − 1)/2, |2 − λ₀| is near 3 and F'_β grows like 3^q
        fits["q*beta"] = fit([r["q"] for r in nodes], [r["q"] * r["beta"] for r in nodes])
    return {"family": name, "rows": rows, "fits": fits}


def main(Q: int) -> None:
    out = [family("p = 1", lambda q: 1, [q for q in range(8, Q + 1, 8)]),
           family("p = (q-1)/2", lambda q: (q - 1) // 2, [q for q in range(9, Q + 1, 8)])]
    OUT.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    for f in out:
        print(f["family"])
        for k, v in f["fits"].items():
            print(f"  {k:14s} a = {v['a']:9.4f}  b = {v['b']:9.2f}  c = {v['c']:10.1f}  max residual {v['max_residual']:.1e}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 256)
