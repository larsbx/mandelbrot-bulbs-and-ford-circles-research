"""The no-limits audit of the finite register: it passes, it bites, and its exemptions are bounded."""

from __future__ import annotations

from pathlib import Path

import pytest

from tools.audit_limits import EXEMPT_SECTIONS, REGISTER, audit

ROOT = Path(__file__).resolve().parents[1]


def document(claims: str, exempt: str = "") -> str:
    """A minimal register: the three exempt sections, and a claims section."""
    first, *rest = EXEMPT_SECTIONS
    parts = [f"# Title\n\n## {first}\n\n{exempt}\n", "## 2. PROVEN\n\n" + claims + "\n"]
    parts += [f"## {heading}\n\nnothing\n" for heading in rest]
    return "\n".join(parts)


def terms(claims: str, exempt: str = "") -> list[str]:
    return [breach.term for breach in audit(document(claims, exempt))]


def test_the_finite_register_passes():
    assert audit(REGISTER.read_text(encoding="utf-8")) == []


def test_the_original_register_would_fail():
    original = (ROOT / "RESEARCH_bulb-ford-correction.md").read_text(encoding="utf-8")
    assert len(audit(document(original))) > 50


@pytest.mark.parametrize(
    ("claim", "term"),
    [
        ("G tends to 1.12 along the family.", "tends to"),
        ("The sequence converges fast.", "converges"),
        ("so `κ → 0.0545` along p = 3.", "→"),
        ("with `lim_{q} G` equal to 1.1.", "lim_"),
        ("d = 2|c'| q⁻² (1 + O(q⁻²)).", "O("),
        ("G = Ĝ(x̃) + o(1).", "o("),
        ("The asymptotic law holds.", "asymptotic"),
        ("The one-sided limit at ⅓.", "limit"),
        ("G approaches 1.008.", "approaches"),
    ],
)
def test_each_limit_idiom_is_caught_in_a_claim(claim, term):
    assert terms(claim) == [term]


def test_convergents_and_convergent_series_are_not_limits():
    assert terms("The convergents p_n/q_n of 3/7; a convergent power series in ℂ{ε}.") == []


def test_a_near_denial_exempts_and_a_far_one_does_not():
    assert terms("No law `d = … (1 + O(q⁻²))` is asserted.") == []
    assert terms("These rows do not support the reading `O(q⁻²)`.") == []
    far = "Nothing here is a problem, and the values at the largest computed N show that G converges."
    assert terms(far) == ["converges"]
    # A denial in a previous clause does not reach into the next one.
    assert terms("No bound is claimed. G converges.") == ["converges"]


def test_the_marker_exempts_exactly_one_paragraph():
    claims = "<!-- limit-exempt: quoting the original -->\nG → 1.12 as quoted.\n\nG → 1.12 again."
    assert terms(claims) == ["→"]


def test_exempt_sections_exempt_and_only_they_do():
    assert terms("A table at q = 59.", exempt="Replaces `lim` and `O(q⁻²)`.") == []


def test_a_missing_exempt_section_refuses_to_run():
    text = document("A table.").replace(f"## {EXEMPT_SECTIONS[-1]}", "## 8. Renamed")
    with pytest.raises(SystemExit):
        audit(text)
