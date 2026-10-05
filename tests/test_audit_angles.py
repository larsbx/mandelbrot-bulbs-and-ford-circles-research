"""The no-angles audit of the kernel: it passes, it bites, and its exemptions are bounded."""

from __future__ import annotations

from pathlib import Path

import pytest

from tools.audit_angles import KERNEL, NUMERICAL_LANE, audit, audit_source


def test_the_kernel_passes():
    assert audit(KERNEL) == []


@pytest.mark.parametrize(
    ("source", "name"),
    [
        ("from math import pi\n", "pi"),
        ("import math\nx = math.pi\n", "pi"),
        ("lam = cmath.exp(2j * pi * p / q)\n", "pi"),
        ("lam = np.exp(2j * np.pi * p / q)\n", "pi"),
        ("lam = acb(arb(2 * p) / q).exp_pi_i()\n", "exp_pi_i"),
        ("w = mp.expjpi(2 * j / n)\n", "expjpi"),
        ("w = mp.expj(t)\n", "expj"),
        ("y = math.sin(t)\n", "sin"),
        ("y = np.cos(t)\n", "cos"),
        ("y = mp.cot(t)\n", "cot"),
        ("y = math.atan2(b, a)\n", "atan2"),
        ("y = math.radians(d)\n", "radians"),
        ("t = arb.pi()\n", "pi"),
    ],
)
def test_angle_idioms_are_caught(source, name):
    assert [b.name for b in audit_source(source, "kernel/bulbford/x.py")] == [name]


def test_prose_and_quadrances_are_not_angles():
    source = '''"""λ₀ = e^{2πip/q}; sin(πp/q) appears here only as prose."""
# Qd(1 − λ₀) = 4 sin²(πp/q), stated in a comment
def quadrance(x, y):
    spin = x * x + y * y        # names that merely contain sin/cos are not calls of them
    cosine_bracket = spin
    return cosine_bracket
'''
    assert audit_source(source, "kernel/bulbford/x.py") == []


def test_exemptions_are_declared_with_reasons_and_exist():
    for path, reason in NUMERICAL_LANE.items():
        assert (KERNEL.parent / path).exists(), path
        assert len(reason) > 20, path


def test_a_stale_exemption_fails(tmp_path: Path):
    pkg = tmp_path / "kernel" / "bulbford"
    pkg.mkdir(parents=True)
    for path in NUMERICAL_LANE:
        (tmp_path / path).write_text("x = 1\n")         # declared numerical, but no angle left
    breaches = audit(tmp_path / "kernel")
    assert {b.name for b in breaches} == {"stale exemption"}
    assert len(breaches) == len(NUMERICAL_LANE)


def test_an_undeclared_module_fails(tmp_path: Path):
    pkg = tmp_path / "kernel" / "bulbford"
    pkg.mkdir(parents=True)
    for path in NUMERICAL_LANE:
        (tmp_path / path).write_text("from math import pi\n")
    (pkg / "exact.py").write_text("from math import pi\nz = pi\n")
    breaches = audit(tmp_path / "kernel")
    assert [(b.path, b.line, b.name) for b in breaches] == [
        ("kernel/bulbford/exact.py", 1, "pi"),
        ("kernel/bulbford/exact.py", 2, "pi"),
    ]
