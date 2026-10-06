<!--
Derived from templates/docs/CONTRIBUTING.md in larsbx/agent-icm @ sha256:88bf9172c22bc8da
Edit the canonical template or estate.toml in larsbx/agent-icm, then re-render there: make estate
Hand-edits here are drift, and agent-icm's `make estate-check` fails on them.
-->

# Contributing to mandelbrot-bulbs-and-ford-circles-research

The size of the satellite bulbs of the Mandelbrot set, and of every hyperbolic
component of the quadratic family, beyond the Ford-circle law q⁻²: a claim
register, its finite (no-limits) version, and finite certificates.

**Language / toolchain:** Python 3.11 (canonical, `kernel/bulbford/`) with numpy, python-flint and
  sympy; pytest
**CI:** GitHub Actions: `ci.yml`, one `test` job: the pinned estate audit, vendored
  digests, the no-limits and no-angles audits, the pytest suite, and the
  vendored digests again

Read these first — they are normative, not background:

- `ESTATE.toml`
- `ARCHITECTURE.md`
- `README.md`
- `RESEARCH_bulb-ford-correction.md`
- `RESEARCH_bulb-ford-correction_finite.md`
- `EXACT_EVIDENCE_BOUNDARY.md`
- `vendored.toml`

---

## The gates

Run these before you open a pull request. Paste what they said into the PR's
evidence table.

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

A check you did not run is not evidence. Say which ones you skipped and why;
the pull request template has a place for exactly that.

## What counts as evidence here

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

These are not style preferences. Each one is settled somewhere in the documents
above; changing one is a decision record, not a pull request comment.

## Working shape

1. **Branch** from the default branch.
2. **Make the failing case first** where this repository's discipline requires
   it, and in every case make sure the new test fails without your change.
3. **Run the gates.** All of them, or name the ones you did not.
4. **Update the surfaces.** Documentation, status tables, ledgers and generated
   artifacts that name the behaviour you changed are part of the change, not a
   follow-up. Regenerate generated files with their tooling; never hand-edit one.
5. **Open the pull request** using the template. Fill in *What this does not
   establish* — it is required, and it is the section reviewers read first.

## Claim discipline

State exactly what your change establishes and no more.

- A search that stopped at a limit reports where it stopped.
- A bounded failure is not an absence.
- A refusal is not a clean answer.
- A translation preserves or lowers authority; it never raises it.
- "Verified" unqualified is not a claim. Say verified *by what*.

## Commits

Imperative, present tense, describing the difference: `Add the M-adic ball
carrier`, `Reject a singular M before the zeroth power`. The body carries the
reasoning when the subject cannot.
