"""P6 sweep: ι_{p/q} = −Σ 1/(1 − F'(z)) over the other fixed points of f^q, for 2 ≤ q ≤ Q.

For every p ≤ q/2 (the rest by conjugation), records the cycle sum against the
Arb ball of `bulbford.index`, the leading cycles, the share of the flank cycle
through α_{w−1} + 1/M and α_{w+2} − 1/M, and the shares of the rotation cycles
and of β. Writes data/cycle_index_sum.json.

    PYTHONPATH=. python3 scripts/cycle_index_sweep.py [Q]
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

from bulbford.cf import coprime_numerators
from bulbford.cycles import cycle_terms, fixed_points, flank_angles, index_by_cycles
from bulbford.index import index_complex

OUT = Path(__file__).resolve().parents[1] / "data" / "cycle_index_sum.json"


def pair(z: complex) -> list[float]:
    return [z.real, z.imag]


def record(p: int, q: int) -> dict:
    terms = sorted(cycle_terms(p, q), key=lambda t: -abs(t.contribution))
    mass = sum(abs(t.contribution) for t in terms)
    flank = set(flank_angles(p, q))
    iota, arb = index_by_cycles(p, q), index_complex(p, q)
    rotation = [t for t in terms if t.rotation is not None and t.period == q]
    return {
        "p": p, "q": q, "pbar": pow(p, -1, q),
        "iota_cycles": pair(iota), "iota_arb": pair(arb), "deviation": abs(iota - arb),
        "cycles": len(terms), "points": 2**q - q - 1, "mass": mass,
        "flank_rank": next(i for i, t in enumerate(terms) if flank <= set(t.angles)),
        "flank_share": abs(terms[0].contribution) / mass,
        "flank_multiplier_modulus": abs(terms[0].multiplier) ** (q // terms[0].period),
        "top10_share": sum(abs(t.contribution) for t in terms[:10]) / mass,
        "rotation_share": sum(abs(t.contribution) for t in rotation) / mass,
        "beta_share": sum(abs(t.contribution) for t in terms if t.period == 1) / mass,
        "top": [
            {"period": t.period, "rotation": None if t.rotation is None else str(t.rotation),
             "min_angle": min(t.angles), "multiplier_modulus": abs(t.multiplier) ** (q // t.period),
             "contribution": pair(t.contribution), "share": abs(t.contribution) / mass}
            for t in terms[:5]
        ],
    }


def main(Q: int) -> None:
    rows = []
    for q in range(2, Q + 1):
        t0 = time.time()
        rows += [record(p, q) for p in coprime_numerators(q) if 2 * p <= q]
        cycle_terms.cache_clear()
        fixed_points.cache_clear()
        print(f"q={q:2d}  pairs={sum(r['q'] == q for r in rows)}  {time.time() - t0:6.1f}s", flush=True)
    OUT.write_text(json.dumps(rows, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {OUT.name}: {len(rows)} p/q, max deviation {max(r['deviation'] for r in rows):.1e}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
