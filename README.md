# mandelbrot-bulbs-and-ford-circles-research

Research branch `bulb-ford-correction`: the size of the satellite bulbs of the Mandelbrot set
(and of every hyperbolic component of the quadratic family) beyond the Ford-circle law `q⁻²`.

- **Register:** [`RESEARCH_bulb-ford-correction.md`](RESEARCH_bulb-ford-correction.md) — PROVEN / VALIDATED / CONJECTURED / FALSIFIED, with the theorem
  `ρ = 1 − q²ε + q³(ι_{p/q} − ½)ε² + O(ε³)` (P2) tying the second-order bulb shape to the holomorphic index of the parabolic point.
- **Finite version:** [`RESEARCH_bulb-ford-correction_finite.md`](RESEARCH_bulb-ford-correction_finite.md) — the same register without limits: exact identities (P0 + P2 as one Laurent-coefficient identity, P6 the projective index sum), finite tables at stated `q`, fitted statistics, Lipschitz-class falsifications, and the conjectures as never-asserted classical referents. The rule is enforced by `tools/audit_limits.py` (CI, `tests/test_audit_limits.py`).
- **Code:** `kernel/bulbford/` (dynamics, parabolic index in ball arithmetic, Taylor data of the multiplier), `tests/` (pytest), `experiments/scripts/` (sweeps), `experiments/data/` (JSON results), `reference/legacy/` (the original instruments, unmodified).
- **Finite certificates:** `kernel/bulbford/certify.py` (Krawczyk + same-box collision exclusions from the finite-Mandelbrot certificate calculus; `experiments/data/center_certificates.json`, 79 satellite centres, `q ≤ 16`), `kernel/bulbford/antipode.py` (certified `ρ = −1` points and `G_ant` enclosures, `experiments/data/antipode_certificates.json`), `kernel/bulbford/wake.py` (exact wake angles, Farey order), and a replay of the `finite-math-kernels` cyclotomic vectors (`tests/vectors/cyclotomic_germ_v1_vectors.json`).
- **Wake visual guide:** [`docs/wake-cycle-to-mandelbrot.md`](docs/wake-cycle-to-mandelbrot.md) shows how the exact doubling cycle selects `θ₋, θ₊`, how the imported rational parameter-ray landing identifies their common bulb root, and where the exact/imported/numerical boundaries lie. The `3/7` still and the interactive browser visual are drawn in [`larsbx/math-vizops`](https://github.com/larsbx/math-vizops) from `kernel/bulbford/wake.py`.
- **Contributions audit:** [`docs/contributions-audit.md`](docs/contributions-audit.md) — what is new, at which tier, against which prior art, and eight consistency findings with the register edits they led to; P2 re-derived independently in `experiments/scripts/audit_p2_independent.py`.
- **Bridges spike:** [`docs/bridges-spike.md`](docs/bridges-spike.md) — four bridges tested on data: the `p̄/q` spectrum (B1), Ramanujan sums and the jump law at every denominator (B2), `2^q − 1` dividing the norm of the parabolic coefficient (B3, `kernel/bulbford/norms.py`), and Dedekind sums (B4, negative).
- **Vendored kernels:** `vendor/python/` holds four packages of [`larsbx/finite-math-kernels`](https://github.com/larsbx/finite-math-kernels), copied byte-for-byte and pinned by SHA-256 in [`vendored.toml`](vendored.toml) (see [Vendoring](#vendoring)): `rational_dynamics_py` (continued fractions, units, Farey sequences, rotation cycles, wakes, rotation numbers, doubling orbits, Dedekind and Ramanujan sums), `vendoring` (the checker), and `lexical_audit` with the `claim_governance` lexer it reads (the engine of `tools/audit_limits.py`, which is this repository's policy over it). `kernel/bulbford/{cf,wake,cycles}.py` are thin adapters over the first.
- **Evidence boundary:** [`EXACT_EVIDENCE_BOUNDARY.md`](EXACT_EVIDENCE_BOUNDARY.md) separates exact low-q coefficient fixtures from ball/FFT/continuation evidence. The canonical low-q corpus is `tests/vectors/parabolic_index_exact_vectors.json`.

```
pip install numpy python-flint sympy pytest
PYTHONPATH=kernel pytest -q
PYTHONPATH=kernel python3 experiments/scripts/analyze.py 59 127 251
PYTHONPATH=kernel python3 experiments/scripts/analyze_dense.py 1009
PYTHONPATH=kernel python3 experiments/scripts/certify_bulbs.py --check
```

## Vendoring

The generic exact arithmetic of `p/q` and of angle doubling is not written here: it is
`rational_dynamics_py`, which `finite-math-kernels` ported from this repository's
`kernel/bulbford/{cf,wake,cycles}.py` and `experiments/scripts/{spectral,bridges_spike}.py`.
`vendor/python/` is a byte-for-byte copy at the commit recorded in `vendored.toml`, and
`vendor/python/vendoring/check_vendored_sync.py` fails on any drift (CI runs it before and
after the tests; `tests/test_vendored.py` runs it too and refuses a local redefinition of a
vendored function). `ESTATE.toml`'s `[[dep]] finite-math-kernels` pin is derived from
`vendored.toml` by the checker, never hand-written. Never edit a vendored file; to update:

```
SHA=<40-hex commit of finite-math-kernels>
rm -rf /tmp/fmk vendor/python/{rational_dynamics_py,vendoring,claim_governance,lexical_audit}   # a clean copy: no stale file survives
mkdir -p /tmp/fmk && git -C ../finite-math-kernels archive $SHA oracles/rational_dynamics_py tools/{vendoring,claim_governance,lexical_audit} | tar -x -C /tmp/fmk
cp -r /tmp/fmk/oracles/rational_dynamics_py /tmp/fmk/tools/{vendoring,claim_governance,lexical_audit} vendor/python/
for name in rational_dynamics_py vendoring claim_governance lexical_audit; do
  python vendor/python/vendoring/check_vendored_sync.py pin $name $SHA   # also re-derives the ESTATE.toml pin
done
python vendor/python/vendoring/check_vendored_sync.py                       # verify
```

`import bulbford` (which `experiments/scripts/paths.py` does) puts `vendor/python` first on
`sys.path`, so `PYTHONPATH=kernel` still suffices, and refuses any other `rational_dynamics_py`
that would shadow the pinned copy. A built wheel carries the vendored package beside `bulbford`.
`wake.py` loaded alone as a file, as `larsbx/math-vizops` does from a sibling checkout, imports
this checkout's `bulbford` first, so the same guard applies.

The adapters keep this repository's contracts where the vendored function differs:
`coprime_numerators(1) == ()` (`units(1) == (0,)`); `farey(n)` is the interior of `F_n`
(`farey_sequence(n, interior=True)`); `modinv` and `xstar` refuse a non-unit `p` instead of
reducing `p/q`. Behaviour that changed with the switch: `farey(0)`, `modinv(p, 1)`, `cf` of a negative `p`,
`from_cf` of an empty or non-positive tail and `mechanical` of a non-reduced `p/q` are now refused; `rotation_number` on a set not closed
under doubling raises `ValueError` (was `KeyError`), and `_orbit` with an even modulus raises
(was an infinite loop); `spectral.ramanujan` is the exact integer `c_q(m)` (was a float cosine
sum, equal to rounding); Dedekind sums use `dedekind_sum`, which divides out `gcd(h, k)` — the
removed `bridges_spike.dedekind` gave `s(2, 4) = −1/32` instead of `0`. Its one caller (B4)
uses the prime `q = 1009`, where the two agree, so `bridges_spike.json` is unchanged.
