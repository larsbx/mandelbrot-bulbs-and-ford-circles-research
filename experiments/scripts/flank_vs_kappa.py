"""Flank cycle against κ and G, on the fractions of experiments/data/cycle_index_sum.json.

For every p/q with 3 ≤ q ≤ 20, p ≤ q/2: the flank |F'| (V22), κ = (ι − ½)/q from
the Arb index, and G_ant from `bulbford.dynamics.bulb`. Reports Spearman rank
correlations within each q (p ≥ 2, where the flank order in x* holds) and
writes experiments/data/flank_vs_kappa.json.

    PYTHONPATH=kernel python3 experiments/scripts/flank_vs_kappa.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from bulbford.cf import xstar
from bulbford.dynamics import MAIN2, bulb
from paths import DATA


def spearman(a: list[float], b: list[float]) -> float:
    rank = lambda v: np.argsort(np.argsort(v)).astype(float)  # noqa: E731
    return float(np.corrcoef(rank(a), rank(b))[0, 1])


def rows() -> list[dict]:
    out = []
    for r in json.loads((DATA / "cycle_index_sum.json").read_text(encoding="utf-8")):
        p, q = r["p"], r["q"]
        if q < 3:
            continue
        kappa = (complex(*r["iota_arb"]) - 0.5) / q
        out.append({"p": p, "q": q, "xstar": float(xstar(p, q)), "flank": r["flank_multiplier_modulus"],
                    "re_kappa": kappa.real, "im_kappa": kappa.imag, "G": bulb(MAIN2, p, q).G_ant})
    return out


def correlations(table: list[dict]) -> list[dict]:
    out = []
    for q in sorted({r["q"] for r in table}):
        s = [r for r in table if r["q"] == q and r["p"] >= 2]
        if len(s) < 3:
            continue
        col = lambda k: [r[k] for r in s]  # noqa: E731
        out.append({"q": q, "n": len(s), **{
            f"{a}~{b}": spearman(col(a), col(b))
            for a, b in [("flank", "re_kappa"), ("flank", "im_kappa"), ("flank", "G"),
                         ("flank", "xstar"), ("re_kappa", "xstar"), ("G", "xstar")]
        }})
    return out


def main() -> None:
    table = rows()
    corr = correlations(table)
    (DATA / "flank_vs_kappa.json").write_text(json.dumps({"rows": table, "spearman": corr}, indent=1) + "\n", encoding="utf-8")
    keys = [k for k in corr[0] if "~" in k]
    print(" q  n  " + "  ".join(f"{k:>16s}" for k in keys))
    for c in corr:
        print(f"{c['q']:2d} {c['n']:2d}  " + "  ".join(f"{c[k]:+16.2f}" for k in keys))


if __name__ == "__main__":
    main()
