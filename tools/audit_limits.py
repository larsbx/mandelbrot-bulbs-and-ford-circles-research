#!/usr/bin/env python3
"""Audit the no-limits rule of the finite register.

`RESEARCH_bulb-ford-correction_finite.md` asserts no limit: every claim is an
exact identity, a certificate, a finite table, a fitted statistic, a finite
falsification, or a classical referent that is recorded and never asserted.
This audit enforces that rule on the prose.

A limit idiom (`→`, `lim`, `O(`, `o(`, "limit", "converges", "asymptotic",
"tends to", "approaches") may appear only where it is covered by one of three
mechanical exemptions:

1. a section whose heading is declared in `EXEMPT_SECTIONS` (the statement
   forms, the classical referents, and the map from the original register:
   the places that name what was removed);
2. a clause that denies the idiom: a denial word precedes it within
   `DENIAL_REACH` characters of the same clause ("No claim below is a limit",
   "do not support the reading `O(q⁻²)`");
3. a paragraph preceded by `<!-- limit-exempt: reason -->`.

"Convergent(s)" is not a limit idiom: continued-fraction convergents and
convergent power series are finite or algebraic objects.
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTER = ROOT / "RESEARCH_bulb-ford-correction_finite.md"

BANNED = re.compile(
    r"→|⟶|->"
    r"|\blim(?:_|\b|sup|inf)"
    r"|\blimits?\b"
    r"|(?<![\w])[Oo]\("
    r"|\bconverg(?:e|es|ed|ence|ing)\b"
    r"|\basymptotic\w*"
    r"|\btends? to\b"
    r"|\bapproach(?:es|ed|ing)?\b",
    re.IGNORECASE,
)

#: Headings (without the leading `## `) of the sections that name what was
#: removed. The audit refuses to run if one is missing, so renaming a section
#: cannot silently widen or narrow the exemption.
EXEMPT_SECTIONS = (
    "0. Statement forms",
    "5. CLASSICAL REFERENTS (never asserted here)",
    "8. Map from the original register",
)

DENIAL = re.compile(r"\b(no|not|nor|never|none|nothing|without|withdrawn|removed|replaces?|unused)\b", re.IGNORECASE)
DENIAL_REACH = 48
CLAUSE = re.compile(r"(?<=[.;])\s+")
EXEMPT_MARKER = re.compile(r"<!--\s*limit-exempt:\s*(?P<reason>[^>]+?)\s*-->")
HEADING = re.compile(r"^(#{1,6})\s+(?P<title>.+?)\s*$")


@dataclass(frozen=True)
class Breach:
    line: int
    term: str
    text: str

    def __str__(self) -> str:
        return f"{REGISTER.name}:{self.line}: limit idiom {self.term!r} outside an exemption: {self.text.strip()[:120]}"


def denied(clause: str, start: int) -> bool:
    return DENIAL.search(clause[max(0, start - DENIAL_REACH):start]) is not None


def clause_breaches(line: str) -> list[str]:
    return [
        match.group(0)
        for clause in CLAUSE.split(line)
        for match in BANNED.finditer(clause)
        if not denied(clause, match.start())
    ]


def audit(text: str) -> list[Breach]:
    lines = text.splitlines()
    headings = {m.group("title") for line in lines if (m := HEADING.match(line)) and len(m.group(1)) == 2}
    missing = [h for h in EXEMPT_SECTIONS if h not in headings]
    if missing:
        raise SystemExit(f"audit_limits: exempt sections missing from {REGISTER.name}: {missing}")

    breaches: list[Breach] = []
    in_exempt_section = False
    marker_armed = False
    in_marked_paragraph = False
    for number, line in enumerate(lines, start=1):
        heading = HEADING.match(line)
        if heading and len(heading.group(1)) == 2:
            in_exempt_section = heading.group("title") in EXEMPT_SECTIONS
            in_marked_paragraph = marker_armed = False
            continue
        if EXEMPT_MARKER.search(line):
            marker_armed = True
            continue
        if not line.strip():
            in_marked_paragraph = False
            continue
        if marker_armed:
            in_marked_paragraph, marker_armed = True, False
        if in_exempt_section or in_marked_paragraph or (heading and number == 1):
            continue
        breaches += [Breach(number, term, line) for term in clause_breaches(line)]
    return breaches


def main() -> int:
    breaches = audit(REGISTER.read_text(encoding="utf-8"))
    for breach in breaches:
        print(breach, file=sys.stderr)
    if breaches:
        return 1
    print(f"no-limits audit passed: {REGISTER.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
