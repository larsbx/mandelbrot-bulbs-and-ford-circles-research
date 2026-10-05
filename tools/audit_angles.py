#!/usr/bin/env python3
"""Audit the no-angles rule of the kernel.

The exact and certified lanes of `kernel/` build every root of unity from its polynomial
(`spread.py`: Sturm on U_{q-1} over Q; `antipode.py`: Krawczyk boxes on X^q - 1;
`index.py`: Arb's roots of Phi_q, selected by order). None of them evaluates an angle,
pi, or a trigonometric function. This audit enforces that on the code itself.

It reads Python tokens, so prose in docstrings and comments (`λ₀ = e^{2πip/q}`) is free;
an identifier token in BANNED is a breach. The modules of the numerical research lane of
EXACT_EVIDENCE_BOUNDARY.md are declared in NUMERICAL_LANE, each with its reason, and are
skipped; a declared module with no angle left is itself a breach, so the exemption list
can only shrink deliberately.

Usage: audit_angles.py    exit 1 and list every breach.
"""

from __future__ import annotations

import io
import sys
import tokenize
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
KERNEL = ROOT / "kernel"

BANNED = frozenset(
    {
        "pi", "tau",
        "sin", "cos", "tan", "cot", "sec", "csc", "sinpi", "cospi", "sinc",
        "asin", "acos", "atan", "atan2", "arcsin", "arccos", "arctan", "arctan2",
        "exp_pi_i", "expj", "expjpi", "cis",
        "radians", "degrees", "deg2rad", "rad2deg",
    }
)

#: Repository-relative paths of the numerical research lane, with the reason each may use angles.
NUMERICAL_LANE = {
    "kernel/bulbford/dynamics.py": "cycle continuation in floating point: a VALIDATED research instrument, seeds only",
    "kernel/bulbford/cycles.py": "floating-point cycle tracking and dense sweeps: VALIDATED numerics, not certificates",
    "kernel/bulbford/taylor.py": "Cauchy/FFT coefficient extraction: the FFT nodes are samples on a circle, VALIDATED",
    "kernel/bulbford/norms.py": "mpmath products over the Galois orbit: numerical cross-checks of exact norms",
    "kernel/bulbford/implosion.py": "Fatou coordinates and horn maps: pi is part of the analytic statements themselves",
}


@dataclass(frozen=True)
class Breach:
    path: str
    line: int
    name: str
    text: str = ""

    def __str__(self) -> str:
        if self.name == "stale exemption":
            return f"{self.path}: declared numerical lane but uses no angle; remove it from NUMERICAL_LANE"
        return f"{self.path}:{self.line}: angle idiom {self.name!r}: {self.text.strip()[:100]}"


def audit_source(source: str, path: str) -> list[Breach]:
    lines = source.splitlines()
    return [
        Breach(path, tok.start[0], tok.string, lines[tok.start[0] - 1])
        for tok in tokenize.generate_tokens(io.StringIO(source).readline)
        if tok.type == tokenize.NAME and tok.string in BANNED
    ]


def audit(kernel: Path = KERNEL) -> list[Breach]:
    root = kernel.parent
    breaches: list[Breach] = []
    for file in sorted(kernel.rglob("*.py")):
        rel = file.relative_to(root).as_posix()
        found = audit_source(file.read_text(encoding="utf-8"), rel)
        if rel not in NUMERICAL_LANE:
            breaches += found
        elif not found:
            breaches.append(Breach(rel, 0, "stale exemption"))
    return breaches


def main() -> int:
    breaches = audit()
    for breach in breaches:
        print(breach, file=sys.stderr)
    if breaches:
        return 1
    print(f"OK: no angle idiom in {KERNEL.relative_to(ROOT)}/ outside the {len(NUMERICAL_LANE)} declared numerical modules.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
