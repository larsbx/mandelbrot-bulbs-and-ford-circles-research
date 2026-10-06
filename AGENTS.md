<!--
Derived from templates/docs/AGENTS.md in larsbx/agent-icm @ sha256:e0ea75600e3d136a
Edit the canonical template or estate.toml in larsbx/agent-icm, then re-render there: make estate
Hand-edits here are drift, and agent-icm's `make estate-check` fails on them.
-->

# Agent policy — mandelbrot-bulbs-and-ford-circles-research

The size of the satellite bulbs of the Mandelbrot set, and of every hyperbolic
component of the quadratic family, beyond the Ford-circle law q⁻²: a claim
register, its finite (no-limits) version, and finite certificates.

**Language / toolchain:** Python 3.11 (canonical, `kernel/bulbford/`) with numpy, python-flint and
  sympy; pytest
**CI:** GitHub Actions: `ci.yml`, one `test` job: the pinned estate audit, vendored
  digests, the no-limits and no-angles audits, the pytest suite, and the
  vendored digests again

This file is for whoever is working here next, human or otherwise. It states
what is settled, so that it does not get re-litigated by someone reading only
the code.

## Read first

- `ESTATE.toml`
- `ARCHITECTURE.md`
- `README.md`
- `RESEARCH_bulb-ford-correction.md`
- `RESEARCH_bulb-ford-correction_finite.md`
- `EXACT_EVIDENCE_BOUNDARY.md`
- `vendored.toml`

## Gates

Before proposing a change as finished, run:

1. vendored packages still match their pins —

   ```sh
   python vendor/python/vendoring/check_vendored_sync.py
   ```

2. the no-limits rule of the finite register —

   ```sh
   python tools/audit_limits.py
   ```

3. the no-angles rule of the kernel —

   ```sh
   python tools/audit_angles.py
   ```

4. suite (CI first installs numpy, python-flint, sympy, pytest, setuptools and
   wheel) —

   ```sh
   python -m pytest -q tests/
   ```

Report honestly which ran. A partial environment that reports a skip is worth
more than one that passes vacuously.

## What this repository treats as evidence

- `RESEARCH_bulb-ford-correction.md` holds claim state: every claim is PROVEN,
  VALIDATED, CONJECTURED or FALSIFIED.
  `RESEARCH_bulb-ford-correction_finite.md` is the same register without
  limits.
- `EXACT_EVIDENCE_BOUNDARY.md` separates exact low-q coefficient fixtures from
  ball, FFT and continuation evidence. The canonical low-q corpus is
  `tests/vectors/parabolic_index_exact_vectors.json`.
- A finite certificate is a box the tests replay from its endpoints alone, in
  rational interval arithmetic with outward dyadic rounding
  (`experiments/data/center_certificates.json`,
  `experiments/data/antipode_certificates.json`). An unmet check is
  INCONCLUSIVE.
- Vendored packages are pinned per file in `vendored.toml`, and the
  `finite-math-kernels` pin in `ESTATE.toml` is derived from it by the
  checker, never hand-written. `tests/test_vendored.py` refuses a local
  redefinition of a vendored function.

## Standing prohibitions

- Never edit a vendored file. Change it upstream in
  `larsbx/finite-math-kernels`, re-vendor and re-pin (README, Vendoring).
- Never assert a limit in the finite register: every claim there is an exact
  identity, a certificate, a finite table, a fitted statistic, a finite
  falsification, or a classical referent recorded and never asserted
  (`tools/audit_limits.py`).
- Never evaluate an angle, pi or a trigonometric function in the exact and
  certified lanes of `kernel/`; roots of unity are built from their
  polynomials (`tools/audit_angles.py`).
- Never read an unmet certificate check as a disproof: it is INCONCLUSIVE, and
  seeds from floating-point continuation are untrusted
  (EXACT_EVIDENCE_BOUNDARY.md).
- Never let a directory rename alone change claim status, acceptance, or
  authority (ARCHITECTURE.md).

## Scope discipline

- Make the change that was asked for. If the surrounding code is wrong in a way
  the task did not name, say so — do not widen the diff to fix it.
- If something is blocked, finish everything that is not, and say precisely what
  was left and why.
- Where a decision is already recorded, follow it or reopen it explicitly. Do
  not route around it in code.
