#!/usr/bin/env python3
"""Audit the no-limits rule of the finite register.

`RESEARCH_bulb-ford-correction_finite.md` asserts no limit: every claim is an
exact identity, a certificate, a finite table, a fitted statistic, a finite
falsification, or a classical referent that is recorded and never asserted.
This audit enforces that rule on the prose. The engine is the vendored
`lexical_audit` package; this file is the policy.

A limit idiom (`→`, `lim`, `O(`, `o(`, "limit", "converges", "asymptotic",
"tends to", "approaches") may appear only where it is covered by one of three
mechanical exemptions:

1. a section whose heading is declared in `EXEMPT_SECTIONS` (the statement
   forms, the classical referents, and the map from the original register:
   the places that name what was removed); a declared section missing from
   the register is itself a finding, so renaming a section cannot silently
   widen or narrow the exemption;
2. a clause that denies the idiom: a denial word precedes it within
   `DENIAL_REACH` characters of the same clause ("No claim below is a limit",
   "do not support the reading `O(q⁻²)`");
3. a paragraph directly preceded by `<!-- limit-exempt: reason -->`.

"Convergent(s)" is not a limit idiom: continued-fraction convergents and
convergent power series are finite or algebraic objects.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "vendor" / "python"))

from lexical_audit import ClauseRule, Policy, Scope, run  # noqa: E402

REGISTER = "RESEARCH_bulb-ford-correction_finite.md"

BANNED = (
    r"(?i)→|⟶|->"
    r"|\blim(?:_|\b|sup|inf)"
    r"|\blimits?\b"
    r"|(?<![\w])[Oo]\("
    r"|\bconverg(?:e|es|ed|ence|ing)\b"
    r"|\basymptotic\w*"
    r"|\btends? to\b"
    r"|\bapproach(?:es|ed|ing)?\b"
)

#: Titles (without the leading `## `) of the sections that name what was removed.
EXEMPT_SECTIONS = (
    "0. Statement forms",
    "5. CLASSICAL REFERENTS (never asserted here)",
    "8. Map from the original register",
)

DENIAL = r"(?i)\b(no|not|nor|never|none|nothing|without|withdrawn|removed|replaces?|unused)\b"
DENIAL_REACH = 48

RULE = ClauseRule(
    scope=Scope((REGISTER,)),
    banned=BANNED,
    message="{path}:{line}: limit idiom {match!r} outside an exemption: {text:.120}",
    denial=DENIAL,
    reach=DENIAL_REACH,
    clause=r"(?<=[.;])\s+",
    exempt_sections=tuple(rf"^{re.escape(title)}$" for title in EXEMPT_SECTIONS),
    marker="limit-exempt",
)

POLICY = Policy(title="no-limits audit", rules=(RULE,))


if __name__ == "__main__":
    sys.exit(run(ROOT, POLICY))
