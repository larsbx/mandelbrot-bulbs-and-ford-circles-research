"""P12: certified continuation from each certified centre (P3) to its certified antipode (P10), q ≤ 16.

Writes data/continuation_certificates.json; `--check` recomputes every path and exits 1 on drift.
Each record is a verdict about the stored P3/P10 boxes: it is recomputed from them, not from floats
(the float Newton only chooses where the next box goes). Runs on all cores (deterministic per bulb).
"""
from __future__ import annotations

import json
import sys
from multiprocessing import Pool
from pathlib import Path

from bulbford.certify import box_from_numerators
from bulbford.continuation import continue_centre_to_antipode

DATA = Path(__file__).resolve().parents[1] / "data"
OUT = DATA / "continuation_certificates.json"


def boxes() -> list[tuple]:
    cen = {(r["p"], r["q"]): r for r in json.loads((DATA / "center_certificates.json").read_text())["certificates"]}
    ant = {(r["p"], r["q"]): r for r in json.loads((DATA / "antipode_certificates.json").read_text())["certificates"]}
    return [(p, q, box_from_numerators(c["box_numerators"], c["prec"]),
             box_from_numerators(ant[(p, q)]["c_box_numerators"], ant[(p, q)]["prec"]),
             box_from_numerators(ant[(p, q)]["z_box_numerators"], ant[(p, q)]["prec"]))
            for (p, q), c in sorted(cen.items(), key=lambda kv: (kv[0][1], kv[0][0]))]


def record(args: tuple) -> dict:
    res = continue_centre_to_antipode(*args)
    return {"p": res.p, "q": res.q, "accepted": res.accepted, "pieces": res.pieces,
            "last_piece_from_s": None if res.last_s0 is None else str(res.last_s0), "reason": res.reason}


def render(workers: int = 4) -> str:
    with Pool(workers) as pool:
        rows = pool.map(record, boxes(), chunksize=1)
    doc = {"schema": "bulbford-continuation-certificates/v1",
           "claim": "the certified antipode lies on the boundary of the hyperbolic component of the certified centre",
           "arithmetic": "Arb complex balls (python-flint acb), exact entry of the rational P3/P10 boxes",
           "certificates": rows}
    return json.dumps(doc, indent=1) + "\n"


if __name__ == "__main__":
    text = render()
    if "--check" in sys.argv[1:]:
        sys.exit(0 if OUT.read_text() == text else 1)
    OUT.write_text(text)
    rows = json.loads(text)["certificates"]
    print(f"wrote {OUT.name}: {sum(r['accepted'] for r in rows)}/{len(rows)} accepted")
