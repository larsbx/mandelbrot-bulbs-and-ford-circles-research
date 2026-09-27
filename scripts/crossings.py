"""The crossings of V25: widest-arc intruders against the ring-2 cycles, q ≤ Q.

For the two families with a flanking arc among the two widest (P9), p = 1 and
p = (q − 1)/2, this computes from single orbits (`bulbford.cycles.cycle_through`)
the terms of the two ring-2 cycles and of every other cycle through an endpoint of
a widest arc (the intruders), and the ratio

    R(q) = max |intruder term| / min |ring-2 term|.

R > 1 means an intruder outranks a ring-2 cycle. Writes data/crossings.json.

    PYTHONPATH=. python3 scripts/crossings.py [Q]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from bulbford.cycles import _orbit, cycle_through
from bulbford.wake import mechanical

OUT = Path(__file__).resolve().parents[1] / "data" / "crossings.json"


def record(p: int, q: int) -> dict:
    M, w, pbar = 2**q - 1, p - 1, pow(p, -1, q)
    A = [mechanical(p, q, r) for r in range(q)]
    flank = set(_orbit((A[(w - 1) % q] + 1) % M, M))
    ring2 = [(A[(w - 2) % q] + 1) % M, (A[(w + 3) % q] - 1) % M]
    taken = flank | {x for a in ring2 for x in _orbit(a, M)} | set(A)
    wide = [x for r in range(q) if (q - 1 if r == q - 1 else (r - w) * pbar % q) >= q - 2
            for x in ((A[r] + 1) % M, (A[(r + 1) % q] - 1) % M)]
    intruders = {}
    for x in wide:
        if x not in taken and not any(x in _orbit(y, M) for y in intruders):
            intruders[x] = abs(cycle_through(p, q, x).contribution)
    r2 = [abs(cycle_through(p, q, a).contribution) for a in ring2]
    return {"p": p, "q": q, "flank": abs(cycle_through(p, q, min(flank)).contribution),
            "ring2": r2, "intruders": sorted(intruders.values(), reverse=True),
            "ratio": max(intruders.values(), default=0.0) / min(r2)}


def main(Q: int) -> None:
    rows = [record(1, q) for q in range(5, Q + 1)] + [record((q - 1) // 2, q) for q in range(7, Q + 1, 2)]
    OUT.write_text(json.dumps(rows, indent=1) + "\n", encoding="utf-8")
    for name, fam in (("p = 1", [r for r in rows if r["p"] == 1]),
                      ("p = (q-1)/2", [r for r in rows if r["p"] == (r["q"] - 1) // 2 and r["p"] > 1])):
        first = next((r["q"] for r in fam if r["ratio"] > 1), None)
        print(f"{name}: first q with R > 1: {first}; R below it: "
              + ", ".join(f"{r['q']}:{r['ratio']:.3f}" for r in fam if first and r["q"] < first)
              + f"; R from it: " + ", ".join(f"{r['q']}:{r['ratio']:.3f}" for r in fam if first and r["q"] >= first)[:400])


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 64)
