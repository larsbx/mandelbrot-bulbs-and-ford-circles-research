# mandelbrot-bulbs-and-ford-circles-research

Research branch `bulb-ford-correction`: the size of the satellite bulbs of the Mandelbrot set
(and of every hyperbolic component of the quadratic family) beyond the Ford-circle law `q⁻²`.

- **Register:** [`RESEARCH_bulb-ford-correction.md`](RESEARCH_bulb-ford-correction.md) — PROVEN / VALIDATED / CONJECTURED / FALSIFIED, with the theorem
  `ρ = 1 − q²ε + q³(ι_{p/q} − ½)ε² + O(ε³)` (P2) tying the second-order bulb shape to the holomorphic index of the parabolic point.
- **Code:** `bulbford/` (dynamics, parabolic index in ball arithmetic, Taylor data of the multiplier), `tests/` (pytest), `scripts/` (sweeps), `data/` (JSON results), `legacy/` (the original instruments, unmodified).

```
pip install numpy python-flint sympy pytest
PYTHONPATH=. pytest -q
PYTHONPATH=. python3 scripts/analyze.py 59 127 251
PYTHONPATH=. python3 scripts/analyze_dense.py 1009
```
