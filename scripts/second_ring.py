"""Second-rank terms of the P6 sum: are they the second ring around the characteristic arc?

With α_0 < … < α_{q−1} the rotation-cycle numerators and w = p − 1 (P7), ring k
consists of the angles α_{w−k} + 1 and α_{w+1+k} − 1 over M = 2^q − 1. Ring 1 is
the flank cycle (P7). For every p ≤ q/2, 5 ≤ q ≤ Q, this records whether the two
ring-2 angles lie on distinct doubling orbits (exact), and whether ranks 2 and 3
of the sum are exactly those two cycles (numerical). For each of the top five
cycles it also records which α-neighbours it contains, as (i − w, ±1).
Writes data/second_ring.json.

    PYTHONPATH=. python3 scripts/second_ring.py [Q]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from bulbford.cf import coprime_numerators
from bulbford.cycles import _orbit, cycle_terms, fixed_points
from bulbford.wake import mechanical

OUT = Path(__file__).resolve().parents[1] / "data" / "second_ring.json"


def ring(p: int, q: int, k: int) -> tuple[int, int]:
    M, A, w = 2**q - 1, [mechanical(p, q, r) for r in range(q)], p - 1
    return (A[(w - k) % q] + 1) % M, (A[(w + 1 + k) % q] - 1) % M


def record(p: int, q: int) -> dict:
    M, A, w = 2**q - 1, [mechanical(p, q, r) for r in range(q)], p - 1
    signed = lambda i: (i - w) % q if (i - w) % q <= q // 2 else (i - w) % q - q  # noqa: E731
    neighbours = {(A[i] + s) % M: [signed(i), s] for i in range(q) for s in (-1, 1)}
    a, b = ring(p, q, 2)
    orbits = frozenset(_orbit(a, M)), frozenset(_orbit(b, M))
    terms = sorted(cycle_terms(p, q), key=lambda t: -abs(t.contribution))
    return {
        "p": p, "q": q,
        "ring2_distinct_orbits": b not in orbits[0],
        "ranks23_are_ring2": {frozenset(terms[1].angles), frozenset(terms[2].angles)} == set(orbits),
        "top": [sorted(neighbours[x] for x in t.angles if x in neighbours) for t in terms[:5]],
    }


def main(Q: int) -> None:
    rows = []
    for q in range(5, Q + 1):
        rows += [record(p, q) for p in coprime_numerators(q) if 2 * p <= q]
        cycle_terms.cache_clear()
        fixed_points.cache_clear()
    OUT.write_text(json.dumps(rows, indent=1) + "\n", encoding="utf-8")
    hits = [r for r in rows if r["ranks23_are_ring2"]]
    print(f"{len(rows)} fractions; ring-2 angles on distinct orbits: {sum(r['ring2_distinct_orbits'] for r in rows)}; "
          f"ranks 2-3 = ring 2: {len(hits)}")
    print("exceptions:", [(r["p"], r["q"]) for r in rows if not r["ranks23_are_ring2"]])


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
