"""The vendored finite-math-kernels packages: pinned, used, and not duplicated.

`rational_dynamics_py` and `vendoring` are copied byte-for-byte into vendor/python and
pinned in vendored.toml. The generic exact arithmetic of p/q and of angle doubling lives
there only: `bulbford.cf`, `bulbford.wake` and `bulbford.cycles` keep their names as thin
adapters that call it, and no other module of kernel/ or experiments/scripts defines one
of its functions again. reference/legacy holds the unmodified original instruments and is
not scanned.
"""
from __future__ import annotations

import ast
import shutil
import subprocess
import sys
import zipfile
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

import rational_dynamics_py as rd
from vendoring import check_vendored_sync as sync

from bulbford import cf, cycles, wake

ROOT = Path(__file__).resolve().parents[1]
SCANNED = sorted((ROOT / "kernel").rglob("*.py")) + sorted((ROOT / "experiments" / "scripts").rglob("*.py"))

#: Names that would duplicate the vendored package: its own, and the local names it replaced.
VENDORED_NAMES = frozenset(rd.__all__) | {
    "modinv", "xstar", "cf", "from_cf", "convergent_denominators", "coprime_numerators",
    "mechanical", "farey", "_orbit", "dedekind", "ramanujan",
}

#: The adapters that keep a local signature, each of which must call the vendored package.
ADAPTERS = {
    "kernel/bulbford/cf.py": {"modinv", "xstar", "cf", "from_cf", "convergent_denominators", "coprime_numerators"},
    "kernel/bulbford/wake.py": {"rotation_cycle", "mechanical", "wake", "farey"},
    "experiments/scripts/spectral.py": {"ramanujan"},
}


def _vendored_aliases(tree: ast.Module) -> set[str]:
    """Names in a module that refer to rational_dynamics_py or one of its functions."""
    out: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module == "rational_dynamics_py":
            out |= {a.asname or a.name for a in node.names}
        elif isinstance(node, ast.Import):
            out |= {a.asname or a.name for a in node.names if a.name == "rational_dynamics_py"}
    return out


def test_vendored_packages_match_their_pins():
    assert sync.repo_root() == ROOT
    assert sync.check() == []
    assert sync.vendored_directories() == ("vendor/python/rational_dynamics_py", "vendor/python/vendoring")


def test_the_vendored_package_is_the_one_imported():
    assert Path(rd.__file__).resolve().parent == ROOT / "vendor" / "python" / "rational_dynamics_py"


def test_no_module_redefines_a_vendored_function():
    for path in SCANNED:
        rel = path.relative_to(ROOT).as_posix()
        tree = ast.parse(path.read_text(encoding="utf-8"))
        aliases = _vendored_aliases(tree)
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) or node.name not in VENDORED_NAMES:
                continue
            assert node.name in ADAPTERS.get(rel, ()), f"{rel}:{node.lineno}: {node.name} duplicates rational_dynamics_py"
            used = {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}
            assert used & aliases, f"{rel}:{node.lineno}: adapter {node.name} does not call rational_dynamics_py"


def test_every_declared_adapter_exists():
    for rel, names in ADAPTERS.items():
        tree = ast.parse((ROOT / rel).read_text(encoding="utf-8"))
        defined = {n.name for n in tree.body if isinstance(n, ast.FunctionDef)}
        assert names <= defined, f"{rel}: stale adapter entries {sorted(names - defined)}"


def test_cycles_reexports_the_vendored_functions():
    assert cycles.rotation_number is rd.rotation_number
    assert cycles._orbit is rd.doubling_orbit


# --- the local contracts the adapters keep, where the vendored function differs ---


def test_coprime_numerators_of_one_is_empty_but_units_of_one_is_zero():
    assert cf.coprime_numerators(1) == () and rd.units(1) == (0,)
    assert all(cf.coprime_numerators(q) == rd.units(q) for q in range(2, 50))


def test_bulbs_farey_is_the_interior_of_the_farey_sequence():
    for n in range(1, 30):
        assert wake.farey(n) == rd.farey_sequence(n, interior=True) == rd.farey_sequence(n)[1:-1]
    with pytest.raises(ValueError):
        wake.farey(0)


def test_modinv_and_xstar_refuse_a_non_unit():
    with pytest.raises(ValueError):
        cf.modinv(2, 4)
    with pytest.raises(ValueError):
        cf.xstar(3, 9)


def test_convergent_denominators_reads_a_non_canonical_expansion():
    assert cf.convergent_denominators((0, 2, 1)) == (1, 2, 3)
    assert cf.convergent_denominators(cf.cf(1, 3)) == (1, 3)


# --- behaviour the vendored package changed ---


def test_dedekind_sum_is_the_defining_sum_when_gcd_exceeds_one():
    """The removed bridges_spike `dedekind` gave s(2, 4) = −1/32; the defining sum is 0."""
    def saw(x: Fraction) -> Fraction:
        return Fraction(0) if x.denominator == 1 else x - (x.numerator // x.denominator) - Fraction(1, 2)

    for k in range(1, 30):
        for h in range(k):
            assert rd.dedekind_sum(h, k) == sum(saw(Fraction(r, k)) * saw(Fraction(h * r, k)) for r in range(1, k))
    assert rd.dedekind_sum(2, 4) == 0


def test_spectral_ramanujan_is_exact():
    from spectral import ramanujan

    m = np.arange(1, 200)
    for q in range(1, 40):
        exact = ramanujan(q, m)
        assert np.array_equal(exact, np.rint(exact))
        assert np.allclose(exact, sum(np.cos(2 * np.pi * a * m / q) for a in rd.units(q)))


def test_rotation_number_and_orbit_refuse_instead_of_failing():
    with pytest.raises(ValueError):
        cycles.rotation_number((1, 2), 7)   # not closed under doubling: was KeyError
    with pytest.raises(ValueError):
        cycles._orbit(1, 8)                 # even modulus: looped forever


def test_adapters_refuse_floats_and_bools_rather_than_coerce():
    for call in (lambda: cf.cf(1.0, 3), lambda: cf.modinv(True, 3), lambda: wake.wake(1, 3.0),
                 lambda: wake.mechanical(1, 3, 0.0), lambda: cycles._orbit(1.0, 7)):
        with pytest.raises(TypeError):
            call()


@pytest.mark.parametrize("q", [1.0, 0.5, True, Fraction(1)])
def test_coprime_numerators_refuses_a_non_integer_even_at_or_below_one(q):
    with pytest.raises(TypeError):
        cf.coprime_numerators(q)


def _import_bulbford(*path: Path) -> subprocess.CompletedProcess:
    """`import bulbford` in a fresh, isolated interpreter whose sys.path starts with `path`."""
    code = (f"import sys; sys.path[:0] = {[str(p) for p in path]!r}; "
            "import bulbford, rational_dynamics_py as rd; print(rd.__file__)")
    return subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True, cwd=ROOT)


@pytest.fixture
def decoy(tmp_path: Path) -> Path:
    """A directory holding another `rational_dynamics_py`."""
    (tmp_path / "rational_dynamics_py").mkdir()
    (tmp_path / "rational_dynamics_py" / "__init__.py").write_text("DECOY = True\n", encoding="utf-8")
    return tmp_path


def test_importing_bulbford_puts_the_pinned_copy_ahead_of_any_other(decoy):
    run = _import_bulbford(decoy, ROOT / "kernel")
    assert run.returncode == 0, run.stderr
    assert Path(run.stdout.strip()).resolve().parent == ROOT / "vendor" / "python" / "rational_dynamics_py"


def test_importing_bulbford_refuses_a_copy_that_shadows_the_pinned_one(decoy):
    run = _import_bulbford(decoy, ROOT / "vendor" / "python", ROOT / "kernel")
    assert run.returncode != 0
    assert "not the pinned copy" in run.stderr


def test_a_built_wheel_carries_the_vendored_package(tmp_path):
    tree = tmp_path / "tree"  # a copy of what the build reads, so the checkout stays clean
    tree.mkdir()
    shutil.copy(ROOT / "pyproject.toml", tree)
    for part in ("kernel", "vendor"):
        shutil.copytree(ROOT / part, tree / part, ignore=shutil.ignore_patterns("__pycache__"))
    build = subprocess.run([sys.executable, "-m", "pip", "wheel", "--no-deps", "--no-build-isolation",
                            "-q", "-w", str(tmp_path), str(tree)], capture_output=True, text=True)
    assert build.returncode == 0, build.stderr
    (wheel,) = tmp_path.glob("bulbford-*.whl")
    names = set(zipfile.ZipFile(wheel).namelist())
    pinned = next(p["files"] for p in sync.load() if p["name"] == "rational_dynamics_py")
    assert set(pinned) <= names
    assert "bulbford/__init__.py" in names


def test_wake_loaded_as_a_bare_file_uses_the_pinned_copy_over_an_ambient_one(decoy):
    code = (f"import sys, importlib.util; sys.path[:0] = [{str(decoy)!r}]; "
            f"spec = importlib.util.spec_from_file_location('bare_wake', {str(ROOT / 'kernel' / 'bulbford' / 'wake.py')!r}); "
            "mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); "
            "import rational_dynamics_py as rd; print(rd.__file__); print(mod.wake(1, 3))")
    run = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True, cwd=decoy)
    assert run.returncode == 0, run.stderr
    found, result = run.stdout.splitlines()
    assert Path(found).resolve().parent == ROOT / "vendor" / "python" / "rational_dynamics_py"
    assert result == "(Fraction(1, 7), Fraction(2, 7))"
