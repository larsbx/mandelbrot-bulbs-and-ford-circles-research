"""Finite certificates for the satellites of the main cardioid, 2 ≤ q ≤ Q_MAX, every unit p.

  data/center_certificates.json    type (0, q) joint boxes for the centres (bulbford.certify)
  data/antipode_certificates.json  two-variable Krawczyk boxes for ρ = −1, with G_ant bounds (bulbford.antipode)

`--check` recomputes both and exits 1 on drift.  Seeds come from the floating-point continuation in
bulbford.dynamics and are untrusted: the tests replay every stored box from its dyadic endpoints alone.
"""
from __future__ import annotations
import cmath, json, sys
from math import pi
from pathlib import Path

from bulbford import antipode, certify
from bulbford.cf import coprime_numerators
from bulbford.dynamics import MAIN2, bulb, center

DATA = Path(__file__).resolve().parents[1] / "data"
Q_MAX = 16
CALCULUS = "finite-mandelbrot-research docs/finite-certificate-calculus.md"
PAIRS = tuple((p, q) for q in range(2, Q_MAX + 1) for p in coprime_numerators(q))


def centre_record(p: int, q: int) -> dict:
    seed = center(MAIN2, MAIN2.c_of(cmath.exp(2j * pi * p / q) * (1 + 1 / q**2)), q)
    return certify.certify_center(p, q, seed).as_record()


def antipode_record(p: int, q: int) -> dict:
    ant = bulb(MAIN2, p, q).ant
    return antipode.certify_antipode(p, q, ant.z, ant.c).as_record()


def document(schema: str, row) -> str:
    doc = {"schema": schema, "calculus": CALCULUS, "family": "z^2 + c, satellites of the main cardioid",
           "tags": certify.TAGS, "certificates": [row(p, q) for p, q in PAIRS]}
    return json.dumps(doc, indent=1, ensure_ascii=False) + "\n"


OUTPUTS = {
    DATA / "center_certificates.json": lambda: document("bulbford-center-certificates/v1", centre_record),
    DATA / "antipode_certificates.json": lambda: document("bulbford-antipode-certificates/v1", antipode_record),
}


if __name__ == "__main__":
    texts = {path: render() for path, render in OUTPUTS.items()}
    if "--check" in sys.argv[1:]:
        sys.exit(0 if all(path.read_text() == text for path, text in texts.items()) else 1)
    for path, text in texts.items():
        path.write_text(text)
        print(f"wrote {path.name}: {text.count('ACCEPTED')} accepted")
