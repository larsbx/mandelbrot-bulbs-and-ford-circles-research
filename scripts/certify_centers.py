"""Joint box certificates (Krawczyk + same-box exclusions) for satellite centres of the main cardioid.

Writes data/center_certificates.json; `--check` recomputes and exits 1 on drift.  Seeds come from the
floating-point continuation in bulbford.dynamics and are untrusted: every stored box is replayed by
tests/test_certify.py from its dyadic endpoints alone.
"""
from __future__ import annotations
import cmath, json, sys
from math import pi
from pathlib import Path

from bulbford.certify import TAGS, certify_center
from bulbford.cf import coprime_numerators
from bulbford.dynamics import MAIN2, center

OUT = Path(__file__).resolve().parents[1] / "data" / "center_certificates.json"
Q_MAX = 16


def seed(p: int, q: int) -> complex:
    return center(MAIN2, MAIN2.c_of(cmath.exp(2j * pi * p / q) * (1 + 1 / q**2)), q)


def render() -> str:
    rows = [certify_center(p, q, seed(p, q)).as_record() for q in range(2, Q_MAX + 1) for p in coprime_numerators(q)]
    doc = {"schema": "bulbford-center-certificates/v1", "calculus": "finite-mandelbrot-research docs/finite-certificate-calculus.md",
           "family": "z^2 + c, satellites of the main cardioid", "tags": TAGS, "certificates": rows}
    return json.dumps(doc, indent=1, ensure_ascii=False) + "\n"


if __name__ == "__main__":
    text = render()
    if "--check" in sys.argv[1:]:
        sys.exit(0 if OUT.read_text() == text else 1)
    OUT.write_text(text)
    print(f"wrote {OUT.name}: {text.count('ACCEPTED')} accepted")
