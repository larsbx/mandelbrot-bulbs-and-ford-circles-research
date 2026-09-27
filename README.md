# mandelbrot-bulbs-and-ford-circles-research

Research branch `bulb-ford-correction`: the size of the satellite bulbs of the Mandelbrot set
(and of every hyperbolic component of the quadratic family) beyond the Ford-circle law `q⁻²`.

- **Register:** [`RESEARCH_bulb-ford-correction.md`](RESEARCH_bulb-ford-correction.md) — PROVEN / VALIDATED / CONJECTURED / FALSIFIED, with the theorem
  `ρ = 1 − q²ε + q³(ι_{p/q} − ½)ε² + O(ε³)` (P2) tying the second-order bulb shape to the holomorphic index of the parabolic point.
- **Finite version:** [`RESEARCH_bulb-ford-correction_finite.md`](RESEARCH_bulb-ford-correction_finite.md) — the same register without limits: exact identities (P0 + P2 as one Laurent-coefficient identity, P6 the projective index sum), finite tables at stated `q`, fitted statistics, Lipschitz-class falsifications, and the conjectures as never-asserted classical referents. The rule is enforced by `tools/audit_limits.py` (CI, `tests/test_audit_limits.py`).
- **Code:** `kernel/bulbford/` (dynamics, parabolic index in ball arithmetic, Taylor data of the multiplier), `tests/` (pytest), `experiments/scripts/` (sweeps), `experiments/data/` (JSON results), `reference/legacy/` (the original instruments, unmodified).
- **Finite certificates:** `kernel/bulbford/certify.py` (Krawczyk + same-box collision exclusions from the finite-Mandelbrot certificate calculus; `experiments/data/center_certificates.json`, 79 satellite centres, `q ≤ 16`), `kernel/bulbford/antipode.py` (certified `ρ = −1` points and `G_ant` enclosures, `experiments/data/antipode_certificates.json`), `kernel/bulbford/wake.py` (exact wake angles, Farey order), and a replay of the `finite-math-kernels` cyclotomic vectors (`tests/vectors/cyclotomic_germ_v1_vectors.json`).
- **Evidence boundary:** [`EXACT_EVIDENCE_BOUNDARY.md`](EXACT_EVIDENCE_BOUNDARY.md) separates exact low-q coefficient fixtures from ball/FFT/continuation evidence. The canonical low-q corpus is `tests/vectors/parabolic_index_exact_vectors.json`.

```
pip install numpy python-flint sympy pytest
PYTHONPATH=kernel pytest -q
PYTHONPATH=kernel python3 experiments/scripts/analyze.py 59 127 251
PYTHONPATH=kernel python3 experiments/scripts/analyze_dense.py 1009
PYTHONPATH=kernel python3 experiments/scripts/certify_bulbs.py --check
```
