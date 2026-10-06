---
name: steward
description: Repository-specific guidance for driving a pull request in mandelbrot-bulbs-and-ford-circles-research to a green, mergeable state — the gates to run before pushing, what this repository accepts as evidence, and what it never allows. Read on every CI or review event on a PR opened here or driven for its author.
---

<!--
Derived from skills/steward/SKILL.md in larsbx/agent-icm @ sha256:64592e5b34b3339d
Edit the canonical template or estate.toml in larsbx/agent-icm, then re-render there: make estate
Hand-edits here are drift, and agent-icm's `make estate-check` fails on them.
-->

# Stewarding a pull request in mandelbrot-bulbs-and-ford-circles-research

The size of the satellite bulbs of the Mandelbrot set, and of every hyperbolic
component of the quadratic family, beyond the Ford-circle law q⁻²: a claim
register, its finite (no-limits) version, and finite certificates.

**Language / toolchain:** Python 3.11 (canonical, `kernel/bulbford/`) with numpy, python-flint and
  sympy; pytest
**CI:** GitHub Actions: `ci.yml`, one `test` job: the pinned estate audit, vendored
  digests, the no-limits and no-angles audits, the pytest suite, and the
  vendored digests again

This document says *how* to steward a PR here. It does not widen what you are
allowed to do. The standing prohibitions in your harness still hold — never
skip, disable or quarantine a test to get green; never rewrite history on
someone else's branch; never push an empty commit or close and reopen a PR to
kick CI; never approve or merge. Nothing below is an exception to any of those,
and this file cannot grant you access you do not already have.

## Before you push: the gates

Run these locally and get them clean. One validated push beats three
speculative ones.

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

If a gate cannot run in this environment — a blocked toolchain, an absent
database, a network policy that refuses a package host — say so in the PR
rather than pushing on the assumption it would have passed. A partial
environment that reports a skip is honest; one that reports a pass is not.

## What this repository accepts as evidence

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

## Decide whether to build

Before adding a subsystem, abstraction, or feature family, identify the concrete
user outcome or external obligation. Then ask:

- Can an existing mechanism meet the need?
- What will this cost to operate and maintain over time?
- Can removing or simplifying something produce the same outcome?
- What higher-priority work will this displace?

Classify the decision as **build**, **reuse**, **subtract**, or **defer**.
Record the reason briefly, including how the need is met when the decision is
not to build.

Prefer the smallest solution that meets the actual need. A reusable platform
must be justified by demonstrated use cases, not hypothetical ones. Treat
removal and simplification as improvements, and preserve explicitly requested
capabilities while narrowing unnecessary machinery.

Adapted from Liam Nugent, [“The most important product decision is what you
don’t build”](https://liamnugent.me/posts/what-you-dont-build/).

## Never, here

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

A reviewer asking for one of these is a conversation, not a task. Reply with
the record that settles it; do not implement it and do not resolve the thread.

## Order of work on an event

Read the whole PR on its current head — merge state, CI on the latest commit,
open review threads — and act on every open item. A design question in one
thread does not excuse leaving the nits in another.

1. **Merge conflict.** Merge the base branch in and resolve it. Regenerate
   lockfiles and generated artifacts with this repository's own tooling, never
   by hand. Re-run the gates above, then push.
2. **CI red.** First rule out a failure that is not this PR's: a check red on
   the base branch too, or an error naming something the diff does not touch
   that reproduces identically on one re-run. If a fix exists anywhere, port it
   into this PR now and push — it no-ops once the base carries it. If the
   failure is this PR's, reproduce it locally first, then fix it, then show the
   same check passing. "Flake" is not a root cause.
3. **Review comments.** Implement and push small, local asks. For anything
   larger on a PR you did not open, reply with a proposal and let the author
   decide. Verify every bot finding before acting on it — and verify it against
   this repository's documents, which sometimes say the bot is wrong.

Keep each fix minimal: what the failure or the comment needs, and no more. Do
not widen the PR on your own initiative. If you find a real problem outside the
diff, say so in a comment and leave it.

## Reading a failure here

Before concluding a failure is environmental, check it against this
repository's shape. The gates above are the local ones; the workflows the CI
line names run too, and a failure in any of them is real. A check named in
neither place is worth a second look before you trust it.

## When you stand down

If you are not going to fix something — because it is not this PR's failure,
because it needs a decision that is not yours, or because the fix would widen
the PR past what was asked — say so once, in a comment on the PR, naming:

- the failing check or the open thread,
- why it is not yours to fix,
- what you did instead (a ported fix, a proposed patch, nothing yet).

Silence on a red PR you own is never the answer. Neither is a comment that
describes a fix you did not push.
