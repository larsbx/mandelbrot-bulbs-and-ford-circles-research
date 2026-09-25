# RESEARCH BRANCH — Bulb geometry of M vs Ford circles: the arithmetic correction Ĝ

Branch id: `bulb-ford-correction`  ·  opened 2026-09-24  ·  revised 2026-09-25  ·  status: **active, one theorem (P2), numerics-validated to q = 2003 (dense) / 16001 (p = 1)**
Discipline: PROVEN / VALIDATED (numerical, stated precision) / CONJECTURED / FALSIFIED are absolute. Analytic imports carry theorem tags. Superseded items are struck, never deleted.

---

## 0. One-paragraph summary

The satellite bulbs of every hyperbolic component of the quadratic family obey
`d(p/q) = 2|φ_H'(e^{2πip/q})| · q⁻² · G(p/q) · (1+O(q⁻²))`, with `G` an `O(1)` correction whose leading dependence is on the modular inverse `x̃ = p̄/q ∈ (−½,½]`, `p·p̄ ≡ 1 (mod q)`: range 0.12 in `|x̃|`, a downward cusp at every rational. **New (P2):** the whole second-order term is controlled by one classical invariant — the holomorphic index `ι_{p/q}` of the parabolic point of `e^{2πip/q}z + z²` — through `ρ = 1 − q²ε + q³(ι − ½)ε² + O(ε³)`; `κ := (ι − ½)/q` is `O(1)`, and `G` is a function of `κ` and the next Taylor coefficient to `10⁻³`. **Correction to the opening statement:** `G` is *not* a function of `x*` alone. A second dependence sits at the *other* end of the continued fraction: bulbs with bounded `p` (`p = 1, 2, 3, …`) lie `+0.018, +0.011, +0.008, …` above the generic `x̃`-limit, generic `t` contributes only `≈ 1e-3`, and this bounded-`p` part is the *only* part that is not universal across hyperbolic components (V19). The generic function `Ĝ(x̃)` runs from `≈1.01` (`x̃ → 0`) to `≈1.16` (near `2−φ`), has cusps in `Re κ` and jumps in `Im κ` at every rational with `q'^{-2}` Ford weights (V21). The Ford-circle analogy is exact at the combinatorial level, asymptotic at leading order, and at first correction the modular group enters through `p ↦ p⁻¹ (mod q)` in the arithmetic of a Galois orbit in `ℚ(ζ_q)`.

---

## 1. Definitions & notation

- `B_{p/q}` : hyperbolic component of period `q` attached to the main cardioid at internal angle `p/q`; root `c_root(p/q) = λ₀/2 − λ₀²/4`, `λ₀ = e^{2πip/q}`. Generic component `H` with multiplier map `λ_H : H → 𝔻`, `c = c_H(λ)`; satellite of `H` at `p/q` has period `q·per(H)`.
- `ε` : `λ = λ₀e^ε`; `u := q²ε` (bulb-scale coordinate). `R_q(u) := ρ(c_H(λ₀e^{u/q²}))`, the satellite multiplier as a function of `u`; `R_q(0)=1`.
- `d(p/q)` : diameter proxy `|c_ant − c_root|`, `c_ant` the boundary point with `ρ = −1` reached from the centre along the root→centre ray. `G(p/q) := d/(2|c_H'(λ₀)| q⁻²)`. `u_a`, `u_c`: solutions of `R_q = −1`, `R_q = 0` near `2`, `1`; `|u_a|/2 = G + O(q⁻²)` (V14).
- `x̃ := p̄/q` reduced to `(−½, ½]`; `x* := |x̃| = q_{n−1}/q = [0; a_n, …, a₁]` for `p/q = [0; a₁, …, a_n]`, `a_n ≥ 2` (tested identity, `tests/test_cf.py`). `t := p/q`.
- `ι_{p/q} := Res_{z=0} dz/(z − f^q(z))`, `f(z) = λ₀z + z²`: holomorphic index of the parabolic fixed point (`= b/a²` for `F = z + az^{q+1} + bz^{2q+1} + …` in equivariant normal form; Écalle's résidu itératif is `(q+1)/2 − ι`). `κ(p/q) := (ι − ½)/q`.
- `r_k(p/q) := [u^k] R_q(u)`; `r₀ = 1`, `r₁ = −1` (P0), `r₂ = κ` (P2).
- Ford circle `C(p/q)`: tangent to ℝ at `p/q`, radius `1/(2q²)`; `C(p/q) ⟂ C(r/s)` tangent iff `|ps−qr| = 1`.

Theorem tags: **[GM84]** Guckenheimer–McGehee `diam B_{p/q} = O(q⁻²)`; **[DH]** Douady–Hubbard multiplier map; **[BS94]** Bullett–Sentenac; **[Yoc95]** Yoccoz; **[MMY]** Marmi–Moussa–Yoccoz; **[Weil]** Kloosterman; **[FL]** Franel–Landau; **[Mil]** Milnor, *Dynamics in One Complex Variable*, §12 (holomorphic index, `ι = b/a²` for `z + az² + bz³`, index sum `Σ 1/(1−ρ) = 0` over finite fixed points of a polynomial); **[BE02]** Buff–Epstein, parabolic Pommerenke–Levin–Yoccoz inequality (bounds `Re resit` for quadratic parabolic points — *statement to be re-checked before use*).

---

## 2. PROVEN (this branch)

- **P0.** First-order local model: `ρ = 1 − q(λ^q − 1) + O((λ^q−1)²) = 1 − q²ε + O(ε²)`; leading law `d = 2|c_H'(λ₀)| q⁻²` to first order. (Elementary.)
- **P1.** `q = 2` exact: period-2 disk diameter `1/2 = 2·sin(π/2)/4`. Also `ρ(λ) = −λ² + 2λ + 4` for the 2-cycle of `λz+z²` (tested).
- **P2 (index formula for the second-order term).** For any unicritical family and any hyperbolic component `H`, with `λ = λ₀e^ε`:
  `ρ(ε) = 1 − q²ε + q³(ι − ½) ε² + O(ε³)`, `ι = ι(f_{λ₀}^{q·per(H)}, z₀)` the holomorphic index of the parabolic cycle point.
  *Proof.* `ρ(λ)` is analytic at `λ₀` (symmetric function of the `q` cycle points, which are the roots of a dynatomic factor inside a fixed disc; residue-theorem continuity). The `q+1` fixed points of `F = f^q` merging at `z₀` — the `H`-cycle point (multiplier `μ = λ^q = e^{qε}`) and the `q` satellite points (multiplier `ρ`) — have total index converging to the index of the merged parabolic point [Mil, 12.9]: `1/(1−μ) + q/(1−ρ) → ι`. With `ρ = 1 − q²ε + Cε² + O(ε³)`: `1/(1−μ) = −1/(qε) + ½ + O(ε)`, `q/(1−ρ) = 1/(qε) + C/q³ + O(ε)`, hence `½ + C/q³ = ι`. ∎
  Checks: `q=1`: `ρ = 2 − e^ε = 1 − ε − ε²/2`, `ι(z+z²) = 0` ✓. `q=2`: `ρ = 1 − 4ε − 3ε²`, `ι(z − 2z³ + z⁴) = 1/8` ✓. Numerically for all `p`, `q ≤ 251` (V12).
- **P3 (integrality, verified exactly for `q ≤ 8`).** With `a = [z^{q+1}] f_{λ₀}^q ∈ ℤ[ζ_q]` (`a = 5 + ζ` at `q=3`, `−10 − 6ζ` at `4`, `36 + 22ζ + 22ζ² + 10ζ³` at `5`, …): `denom(ι) = N(a)²` for `q = 3` and `denom(ι·a²) = 3, 2, 5, 3, 7, 1` for `q = 3…8` — exactly the odd part of `q` (`data/exact_denominators.txt`, `scripts/exact_index.py`). So `q·ι·a² ∈ ℤ[ζ_q]`. *Sketch for general `q` (C22):* `qι = b₂/b₁² + (q−1)/2` with `b₁ = a` (P2 remark (i)); the conjugacy to the equivariant normal form has denominators `λ₀^{j} − λ₀ = λ₀(ζ^{j−1} − 1)`, cyclotomic units unless `q` is a prime power `ℓ^k`, where they generate the prime above `ℓ`; hence `b₂ ∈ ℤ[ζ_q][1/ℓ]` and the `ℓ`-power is bounded by `q`.
  *Remarks.* (i) The `ε`-derivatives of the normal-form coefficients cancel at second order (checked by direct normal-form expansion: `ρ = 1 − q²ε + q²(b₂/b₁² − ½)ε² + …` for `f ~ λw(1 + b₁w^q + b₂w^{2q})`, and `b₂/b₁² = qι − (q−1)/2`). (ii) `ι_{p/q} ∈ ℚ(ζ_q)` and `ι_{p'/q} = σ(ι_{p/q})` for the Galois automorphism `ζ ↦ ζ^{p'p̄}` (the residue is a rational function of `λ₀` with rational coefficients). (iii) Third and higher orders depend on the unfolding, not on the germ alone.

---

## 3. VALIDATED (numerical; doubles unless stated; Newton to 1e-14 relative)

Legacy items V1–V11 stand (reproduced 2026-09-25 by the new instruments: V5/V9/V10 table at `q=59` identical to 4 digits). Corrections: V6's "one-sided limits agree" is superseded by V16; V4's "residual dependence on `a_{n−1}` ≈ 0.005" is explained by V17–V18.

- **V12. P2 verified.** `[u²]R_q` (Cauchy/FFT on `|u| = 2.5`, 256 points) vs `(ι−½)/q` (ball arithmetic, radius `< 1e-12`): max deviation over all `p`: `9e-15 (q=59)`, `1e-13 (127)`, `1.3e-12 (251)`. `r₀ = 1`, `r₁ = −1` to `1e-11`.
- **V13. Exact values** (`q ≤ 12`, `data/exact_index.txt`): `ι_{1/3} = (92 − 16ζ₃)/441`, `ι_{1/4} = (1447 − 365i)/4624`, `ι_{1/5} = (1068·31 − 9153ζ − 7113ζ² − 312ζ³)/93775`; denominators are squares of algebraic norms (`441 = 21²`, `4624 = 68²`, `3025 = 55²`), as expected for `b/a²`. `Im ι(1/q)` decreases by `0.064 ± 0.001` per unit `q` for `q = 3…12`: linear growth, i.e. `κ = O(1)` already visible exactly.
- **V14. The limit shape lives in the `ε`-plane.** `|u_a|/2` (root of the truncated Taylor series of `R_q` at `u=0`, radius of the sampling circle 2.5) reproduces `G_ant` (direct antipode continuation) to `5.7e-4 (q=59)`, `1.3e-4 (127)`, `5.5e-5 (251)` — i.e. `O(q⁻²)`, the nonlinearity of `c_H(λ)`. `|r_k|` decays like `q^{-k}`-ish (tail `< 1e-16` at `k=16`); `R_q` is analytic on `|u| ≲ cq`.
- **V15. `κ` is `O(1)`, and `G` is a function of `(κ, r₃)`.** All `p`, `q ∈ {59,127,251,1009}`: `Re κ ∈ [0.015, 0.064]`, `|Im κ| ≤ 0.055`, `|r₃| ≤ 2.3e-3`, `|r₄| ≤ 1e-4`. `Re κ` is a smooth unimodal function of `x*` with the profile of `Ĝ` (`0.024` at `x*→0`, max `≈0.062` near `0.40`, `0.050` at `½`). `Im κ` is **odd in the signed inverse** `x̃` (`κ(−x̃) = conj κ(x̃)`, exact by conjugation), `|Im κ| ≈ 0.053` at `x*→0`, crossing `0` near `x* ≈ 0.32–0.39`, `0.024` at `½`. Newton on `1 − u + κu² + r₃u³ = −1` from `u=2` gives `G` to `≤ 3e-3` (peak `1.147` ✓, `x*→0`: `1.029` vs `1.026`).
- **V16. One-sided limits at rationals.** Two kinds of sequence give two kinds of limit.
  *(a) Bounded `p`* (`x̃ = [0; tail, N]` exactly, `q ≤ 1153`, converged to 1e-5):
  `⅓⁻` via `[0;3,N]` (`p=3`): `κ → 0.05453 + 0.00506i`, `G → 1.12837`; `⅓⁺` via `[0;2,1,N]` (`p=3`): `κ → 0.05504 − 0.02351i`, `G → 1.12537`;
  `½⁻` via `[0;2,N]` (`p=2`): `κ → 0.04969 + 0.0240i`, `G → 1.11142`; `0⁺` via `[0;N]` (`p=1`): `κ → 0.02383 − 0.0527i` (Im log-slow), `G → 1.0263–1.0265`.
  *(b) Generic* (`x̃ = [0; tail, N, a₁]`, `a₁ ∈ {16, 64}` — the `a₁`-dependence is `< 3e-4` — `N ∈ {64,128,256}`, `q ≤ 16385`; approach is `c/N`, `c ≈ 0.2–0.6`, quoted with the fit; `data/generic_limits.json`):
  `⅓⁻`: `G = 1.12345, 1.12196, 1.12122 → 1.1205`; `κ → 0.0517 + 0.0012i`.
  `⅓⁺`: `G = 1.12152, 1.12002, 1.11927 → 1.1185`; `κ → 0.0521 − 0.0195i`.
  `½⁻`: `G = 1.10505, 1.10257, 1.10133 → 1.100`; `κ → 0.0447 + 0.0173i`.
  `0⁺`: `G = 1.01756, 1.01293, 1.01067` (and `1.00998` at `q=16385`) `→ ≈1.008 (±0.002)`; `κ → ≈0.0116 − 0.0380i`.
  Reading: the generic limits are the legacy V6 values (`1.120/1.121` at `⅓`, `≈1.102` at `½`); the bounded-`p` sequences sit `+0.008 (⅓⁻)`, `+0.007 (⅓⁺)`, `+0.011 (½⁻)`, `+0.018 (0⁺)` above them. At `⅓` the generic `Im κ` **jumps** by `−0.021` (bounded-`p`: `−0.029`), `Re κ` is continuous to `5e-4`, and `G` has a cusp plus a jump of `≈ −0.002 ± 0.001`.
- **V17. The limit depends on the far end of the continued fraction (order of limits).** Same `x̃` to `1e-5`, different `G`:
  `x̃ = [0;3,N,5]` (`p/q = [0;5,N,3]`, `t→⅕`): `G = 1.12551, 1.12462, 1.12371, 1.12326` at `N = 48, 64, 96, 128`, fit `G_∞ + 0.17/N` ⇒ `G_∞ ≈ 1.1219`, vs `1.12837` for `[0;3,N]`.
  `x̃ = [0;N,4]` (`p/q=[0;4,N]`, `t→¼`): `G → ≈1.014` (`1.01509` at `q=1537`, slope `0.5/N`) vs `1.0265` for `p=1`.
  At fixed `N`, `a₁ → ∞` (`p` fixed `= 3N+1`, `q → ∞`) converges smoothly: `[0;3,16,a₁]`: `G = 1.13089, 1.13164, 1.13179, …, 1.13143` for `a₁ = 2 … 100`; `[0;64,a₁]`: `1.02585 → 1.01689`. Grid `[0;3,N,a₁]` (`data/limit_grid.json`): `lim_{a₁→∞} G = 1.12614, 1.12322, 1.12173` for `N = 32, 64, 128` (`c/N`, `c = 0.19`) `→ 1.1202` as `N → ∞`, against `1.12837` for `a₁ = ∞` taken *first* (`p=3`). The `a₁`-dependence at fixed `N` is monotone (`a₁ = 2`: `1.12855 → 1.12614` at `N=32`) and of size `0.002–0.005`; the two orders of limits differ by `0.008`.
- **V18. The second variable is arithmetic in `t`, not smooth.** At tail `(N=128, 3)` (`x̃ = ⅓ − 9e-4`), prefixes giving `t = 0.5, 0.4, 0.4286, 0.4545, 0.333, 0.286, 0.25, 0.2, 0.111`: `G = 1.1260, 1.1238, 1.1232, 1.1241, 1.1246, 1.1236, 1.1238, 1.1233, 1.1238` — peaks at `t = ½` (`+0.0024`), `⅓` (`+0.001`) over a `≈1.1236` floor, and `t→0` with `p=3`: `1.1284`. Symmetry `t ↔ 1−t` exact (`(2,2)≡(1,1,2)`, `(3,)≡(1,2)`, `(4,)≡(1,3)`).
  Dense `q = 1009` (all 504 `p ≤ q/2`, `data/kappa_q1009.json`): `var G = 1.64e-3`; a 100-bin function of `x*` alone leaves residual variance `2.65e-6` (`R² = 0.9984`, residual sd `1.6e-3`); adding a 50-bin function of `t` (backfitting) leaves `1.87e-6` (`R² = 0.9989`, `sd H = 1.0e-3`). At `q = 2003` (1001 bulbs): `R² = 0.9985` (`x*` only), `0.9986` with `t`, `sd H = 4.5e-4` — the generic `t`-part **halves** from `q = 1009` to `2003`, so it is a finite-`q` effect vanishing in the limit; only the bounded-`p` part survives (`p = 1,2,3`: `+0.0081, +0.0059, +0.0021` at `q = 2003`, unchanged). The large residuals (local: `G` minus its 8 nearest `x*`-neighbours) are the bounded-`p` bulbs: `p = 1,2,3,4`: `+0.0080, +0.0044, +0.0025, +0.0005`; and the bulbs with `p·p̄ = q−1` (`x̃ = −1/p` exactly, `6 ≤ p ≤ 36`): `−0.0036 … −0.0053` (they sit at the bottom of the `1/p` cusps). Small-denominator `t` classes (`t ≈ ⅓, ¼, ⅖`, `n = 4` each) average `−0.002 … −0.003` (2–3σ). Generic `p ≥ 20` off these sets: `|e| < 1e-3`.
- **V21. `Im κ` jumps at rationals `p'/q'` (q = 1009, means over `0.002 < |x̃ − p'/q'| < 0.009` on each side, `J := Im κ(p'/q'⁺) − Im κ(p'/q'⁻)`):**
  `½: −0.034 (from the two conjugate sides)`, `⅓: −0.0186`, `¼: −0.0105`, `⅕: −0.0045`, `⅖: −0.0053`, `⅙: −0.0029`, `⅐: −0.0004`, `2/7: −0.0007`, `3/7: 0.0000`.
  `J·q'² = −0.137, −0.167, −0.167, −0.112, −0.132, −0.106` for `q' = 2,3,4,5,5,6`: Ford-type `q'^{-2}` weight, equal for `⅕` and `⅖`, then falling faster at `q' = 7` (window-limited at this `q`). Same sign everywhere (the jump is *down* in `Im κ` for `x̃ > 0`; odd in `x̃`). Together with V16(b): `Re κ` has the cusps, `Im κ` the jumps. *Refined by V23.*
- **V22. `ι` is a collective sum over the repelling `q`-cycles, not dominated by one cycle** (`scripts/cycle_decomposition.py`, `data/cycle_decomposition.txt`). Index sum [Mil]: `ι = −1/(1−(2−λ₀)^q) − q Σ_C 1/(1−ρ_C)` over the other `(2^q−2−q)/q` cycles of period `q`. All `2^q` fixed points of `f^q` found (Newton, 80 000 starts; counts `17/17` at `q=7`, `185/185` at `1/11`), numeric sum `= ι` (series) to 6 digits. Largest single term (least repelling cycle): `1/7`: `|ρ| = 32.4`, term `0.152 − 0.160i` of `ι = 0.508 − 0.261i`; `1/11`: `|ρ| = 41.4`, term `0.179 − 0.202i` of `0.682 − 0.521i`, with the remaining 178 cycles summing to `0.519 − 0.082i`. Consequence for C17′: the asymptotics of `ι` is a partition-function-type statement over `~2^q/q` cycles (the Farey-neighbour picture in next-move 1 is *not* literal for prime `q`, where no cycle of period `q' ≠ q` is a fixed point of `f^q`).
- **V23. Cusp and jump laws at `q = 1009`** (`scripts/cusp_fit.py`, window `0.0025 < |x̃ − p'/q'| < 0.012`, two-sided linear fits, bounded-`p ≤ 5` excluded):
  `G` cusp depth (envelope at `0.012` minus mean intercept) `×q'²`: `⅓: 0.157`, `¼: 0.210`, `⅕: 0.196`, `⅖: 0.245`, `⅙: 0.133`, `⅐: 0.042` — the `≈0.18–0.25 / q'²` law of V7 holds for `q' ≤ 5`, resolution-limited beyond. V-slopes of `G`: `±1.45 (⅓)`, `0.78/1.41 (¼)`, `0.45/0.86 (⅕)`. `Re κ` cusp depth `×q'²`: `0.053, 0.073, 0.070, 0.080, 0.049` for `⅓, ¼, ⅕, ⅖, ⅙`.
  `Im κ` jump from fitted intercepts: `⅓: −0.0216`, `¼: −0.0148`, `⅕: −0.0107`, `⅖: −0.0121`, `⅙: −0.0068`, `⅐: −0.0045`, `2/7: −0.0052`, `3/7: −0.0054`, `3/8: −0.0057`; `J·q' = −0.065, −0.059, −0.054, −0.060, −0.041, −0.032, −0.037, −0.038, −0.046`. With intercept extrapolation the jump decays like `q'^{-1}` (slowly varying prefactor), *not* `q'^{-2}` as the cruder window means of V21 suggested. **Settled at `q = 2003`** (`data/cusp_fit_2003.txt`, window `0.008`, 1001 bulbs): `J = −0.0210 (⅓), −0.0146 (¼), −0.0104 (⅕), −0.0120 (⅖), −0.0092 (⅙), −0.0046 (⅐), −0.0068 (2/7), −0.0067 (3/7), −0.0053 (3/8)` — within `1e-3` of the `q = 1009` values, i.e. converged; `J·q' = −0.063, −0.058, −0.052, −0.060, −0.055, −0.032, −0.048, −0.047, −0.042`. **Law: `J(p'/q') ≈ −(0.05 ± 0.01)/q'`**, depending on `q'` only (`⅕ ≈ ⅖`, `2/7 ≈ 3/7`). `G` cusp depths at window `0.008`, `×q'²`: `0.154 (⅓), 0.203 (¼), 0.159 (⅕), 0.263 (⅖), 0.114 (⅙), 0.098 (⅐), 0.073 (2/7), 0.209 (3/7)` — the `q'^{-2}` Ford weight of V7 holds through `q' = 7` once resolved. The V-slopes steepen as the window shrinks (`⅓`: `1.45 → 2.14` from `0.012` to `0.008`): the cusp is **not linear**. Profile fits `y0 + s|x−x0|^α` (`scripts/cusp_shape.py`, window `0.03`) are contaminated by sub-cusps except at `½`, where `α = 0.80 (G), 0.76 (Re κ)` at both `q = 1009` and `2003`; at `⅓⁺` `α ≈ 0.34–0.42`, elsewhere the fit pins to the grid edge — open (next move 2). `G` intercept jumps stay `≤ 0.005` with no consistent sign: the generic jump of `G` at rationals is `≲ 0.002`, i.e. `Ĝ` is (numerically) continuous while `Im κ̂` is not.
- **V24. The `p = 1` constant (C20).** `κ(1/q)` for `q = 1009, 2003, 4001, 8009, 16001`: `Re = 0.0238263, 0.0238260, 0.0238259, 0.0238259, 0.0238260` (converged); `Im = −0.0524595, −0.0523834, −0.0523443, −0.0523246, −0.0523150`, differences halving per doubling (`c/q`, `c = 0.155`); Richardson `2κ(2q) − κ(q)`: `−0.0523049, −0.0523054` ⇒
  `κ₀ = 0.0238260 − 0.052305i (±1e-6)`. Not `1/42 − …` (`1/42 = 0.0238095`, off by `1.6e-5`). `r₃(1/q) → −0.000377 − 0.000979i`, `G(1/q) → 1.026411`.
- **V26. P2 holds in every family tested.** `index_at(fam, p, q)` (index of the parabolic cycle point of the satellite root, ball arithmetic) vs `[u²]R_q`: `|Δ| < 1e-7` for `disk2` at `2/7, 3/11` (period `2q` cycle) and `main3` (`z³+c`) at `2/7, 3/11` (`tests/test_index.py::test_P2_other_families`).
- **V27. Degree 3 (`z³+c` main component), same sequences as V16** (`scripts/main3_limits.py`): bounded-`p` / generic limits `G₃`: `⅓⁻`: `1.2153 / 1.198`; `½⁻`: `1.2396 / 1.2347`; `0⁺`: `1.175 / ≈1.13` (still falling at `q = 2049`). `κ₃`: `⅓⁻ generic 0.1488 − 0.0014i`, `½⁻ generic 0.155 + 0.121i`, `0⁺ generic 0.110 − 0.027i`. Ratio `G₃/G₂` runs `1.067 (⅓⁻)`, `1.120 (½⁻)`, `1.126 (0⁺ generic)`, `1.145 (0⁺, p=1)`: `Ĝ₃` is a different function (V10 confirmed at the limit level), `Ĝ₃(0⁺) ≈ 1.13 ≠ 1` whereas `Ĝ₂(0⁺) ≈ 1.01`; the bounded-`p` excess is larger in degree 3 (`+0.017` at `⅓⁻` vs `+0.008`).
- **V28. Hölder exponent ≈ ½ (`q = 1009` and `2003`).** Modulus of continuity on the 499 generic bulbs (`scripts/holder.py`), `ω(δ) = max |y(x) − y(x')|`, `|x − x'| ≤ δ`, `δ ∈ [0.003, 0.048]`: `ω_G ∝ δ^{0.59}`, `ω_{Re κ} ∝ δ^{0.52}` at `q = 1009`; `δ^{0.63}`, `δ^{0.45}` at `q = 2003` (996 generic bulbs); `ω_{Im κ}` saturates at the jump sizes (`0.083` beyond `δ = 0.003` at `q = 2003`: the `½` jump). Below `δ ≈ 0.003` the finite-`q` floor (`≈ 0.008` in `G`, unchanged between the two `q`) dominates — that floor is the cusp structure below the sampling scale, not noise. Reading: `Ĝ` and `Re κ̂` have Hölder exponent `≈ 0.5–0.6` — the [MMY] exponent — and are not Lipschitz; consistent with V23's non-linear cusps.
- **V29. `log|a|/q` is Brjuno-like** (`scripts/brjuno_test.py`, `q = 59, 127, 251`, all `p ≤ q/2`; `log|a|` in ball arithmetic). With `B_fin(x) := Σ_{k<n} β_{k−1} log(1/α_k)` (Brjuno sum of the rational `x`, infinite last term dropped): `corr(log|a|/q, B_fin(p/q)) = 0.932, 0.930, 0.932`; fit `log|a|/q = c_q + s_q·B_fin(p/q)` with `(c_q, s_q) = (0.343, 0.418), (0.274, 0.470), (0.221, 0.505)` (`s_q` still rising, `≈ 0.5–0.6` extrapolated); rms residual `0.096, 0.123, 0.138` vs spread `0.266, 0.336, 0.382`; adding `B_fin(x̃)` lowers the residual by `10–20 %`. `corr` with `x*` alone is `−0.25, −0.18, −0.13`. So the *cycle radius* `e^{−log|a|/q}` is governed by the whole continued fraction of `p/q` (Yoccoz-type: `log r = −B + O(1)`), while the *shape correction* `κ` is governed by its last digits (`x̃`) — two different arithmetic functionals of the same fraction.
- **V31. `⟨G⟩ = 1.1122`** (mean over all `p ≤ q/2`: `1.11224` at `q = 1009`, `1.11221` at `2003`; generic-only `1.11236 / 1.11224`; sd `0.0405`). Since `x̃` is equidistributed [Weil], this is `∫Ĝ` and fixes the constants of C2′/C9′: `Res_{s=1} Z_M = (12/π³)⟨G⟩ = 0.4304`, `Σ_{q≤Q} d ≈ (24/π³)·1.1122·log Q = 0.861·log Q` — inside the legacy V11 slope range `0.83–0.89` (which had estimated `⟨Ĝ⟩ ≈ 1.08–1.09` from `Q ≤ 60`, where `G` is still below its limit).
- **V30. Along a canonical family the cusp is a straight line; the Hölder-½ modulus is an envelope over families.** `scripts/canonical_profile.py` / `canonical_fit.py` (`data/canonical_fit.txt`): `x̃ = [0; tail, N, a]`, `N = 16…384`, `q` up to `18451`, `δ = |x̃ − p'/q'| ∈ [3e-4, 8e-3]`. Linear fits `G = G_tip + s·δ` have rms `≤ 1e-4` (`≤ 1e-5` at `½`, `0`):
  `⅓⁻`: `(G_tip, s) = (1.12056, 1.65)` for `a=16`, `(1.12199, 1.49)` for `a=5`; `⅓⁺`: `(1.11857, 1.72)` (`a=16`); `½⁻`: `(1.10008, 1.29)`, `(1.10230, 1.19)`; `0⁺`: `(1.00842, 0.58)`, `(1.01228, 0.54)`. `Re κ` is linear too (`κ_tip` `0.0517 (⅓⁻)`, `0.0522 (⅓⁺)`, `0.0447 (½⁻)`, `0.0116 (0⁺)`; slopes `0.58, 0.58, 0.46, 0.28`); `Im κ` is flat at its one-sided value (`+0.0012 / −0.0195` at `⅓∓`, `+0.0174` at `½⁻`, `−0.038` at `0⁺`).
  Readings: (i) one-sided derivatives of `Ĝ` and `Re κ̂` exist along every family, so the cusps are V-shaped *within* a family; the steepening and the `δ^{0.5–0.6}` modulus of V23/V28 come from the spread of the tips across families (the deep digit `a` shifts the tip by `0.0014 (⅓⁻)`, `0.0022 (½⁻)`, `0.0039 (0⁺)` between `a = 5` and `16`, `< 3e-4` between `16` and `64`); (ii) the generic jump of `Ĝ` at `⅓` is `G_tip(⅓⁺) − G_tip(⅓⁻) = −0.0020` (`a = 16`, rms `1e-4`): small but real; (iii) `Ĝ(0⁺) ≈ 1.008`, `Ĝ(½⁻) ≈ 1.100`, `Ĝ(⅓∓) ≈ 1.1206 / 1.1186` (`a → ∞` extrapolated to `±5e-4`). This supersedes the `c/N` extrapolations of V16(b) (same numbers, now with a fitted linear law).
- **V25. Raw leading coefficients are not the invariant.** `a = [z^{q+1}] f^q`, `b = [z^{2q+1}] f^q` in the coordinate `f = λ₀z + z²` (`scripts/leading_coeff.py`, `q = 59`): `|b/a²|/q` ranges over `10^{10}…10^{30}` while `ι = O(q)` — the non-resonant intermediate terms carry almost all of `b`, so any small-denominator analysis must be done in the equivariant normal form (next move 1), not on raw Taylor data. Clean by-product: `log|a|/q ≈ 1.05–1.16` for generic `p` (rising to `1.27` at `x* = 0.05`), and `2.32, 1.78, 1.52, 1.37, 1.28` for `p = 1…5`; since the satellite cycle radius is `≈ (q|ε|/|a|)^{1/q} → e^{−log|a|/q}`, the cycle sits at radius `≈0.33` generically and `≈0.10` for `p = 1`. Bounded `p` is again the exceptional set.
- **V19. Universality across components is exact for the `x̃`-part and fails only for bounded `p`.** `disk2/main2` ratio of `G` on `[0;a₁,16,3]`: `0.9924` at `p=3` (`a₁=1`), `1.0000 ± 0.0002` for `a₁ = 2…100`; on `[0;a₁,64]`: `0.9881` at `p=1`, `0.9997–1.0008` for `a₁ ≥ 2`. Refines V9 (`0.9985 ± 0.0032` at `q=59`, whose scatter was the bounded-`p` bulbs).
- **V20. Quadratic irrational:** `x̃ → 2−φ = [0;2,1,1,…]`: `κ = 0.05716, 0.06216, 0.06390` at `q = 34, 89, 610`; geometric convergence (ratio `≈ 1/3` per two Fibonacci steps), `κ_∞ ≈ 0.0648 − 0.0058i`, `G_∞ ≈ 1.158` (the global maximum region).

Precision caveats: doubles for `R_q`; the FFT route is accurate to `~1e-12` in `r₂` (verified against balls, V12) and is the instrument for `q > 300`. `G` from `|u_a|/2` carries the `O(q⁻²)` proxy difference of V14 (irrelevant for limits). Prime `q` used in sweeps; sequences use whatever `q` the CF gives.

---

## 4. FALSIFIED (this branch)

- F1 (was C12). `c₁(t)/q` correction with smooth `c₁`. **False** (V2–V3).
- F2 (was C13). "Nothing arithmetic transfers." **False** (V3–V6, V15).
- F3 (was C15). Same `Ĝ` for every degree. **False** (V10).
- **F4 (was C1″/C14′ as stated).** "`d = 2|φ'|q⁻² Ĝ(x*)(1+o(1))` with a single continuous function `Ĝ` on `(0,½]`." **False at the 0.5–1.8 % level**: the `q→∞` limit along `x̃ → x` depends on whether `p` stays bounded (V16–V17: `1.0265` vs `≈1.008` at `0⁺`, `1.1284` vs `1.1205` at `⅓⁻`), and one-sided generic limits at `⅓` differ by `≈0.002` (jump on top of the cusp). The leading `x*`-description survives as the dominant term (range `0.15` vs generic residual `1.6e-3`, R² = 0.9984 at `q = 1009`).
- **F5 (was the contraction form of C16).** "Deeper ancestors enter with weight `(q_{n−1}/q)²`." **False**: the digit `a₁` of `p/q` (deepest ancestor of `x̃`) changes the limit by `0.005–0.012` even as the intermediate digit `N → ∞` (V17), which would be suppressed by `(1/N)²` under C16.
- ~~V6's claim of equal one-sided limits at ⅓~~ (superseded; the `0.001` agreement there was at one fixed prefix and inside its error).

---

## 5. CONJECTURED (open register)

**C1‴ (size law, corrected form).** For every hyperbolic component `H` of the quadratic family,
`d_H(p/q) = 2|c_H'(λ₀)| q⁻² · G_H(p/q) · (1 + O(q⁻²))`, with
`G_H(p/q) = Ĝ(x̃) + Ĥ_H(p; x̃) + o(1)`, where `Ĝ` is universal (component-independent), depends on `x̃ = p̄/q` only, range `≈[1.01, 1.16]`, and is the limit along any sequence with `p → ∞`; `Ĥ_H` is supported on bounded `p` (`+0.018, +0.011, +0.008, +0.0005` for `p = 1,2,3,4` at `x̃ → 0, ½, ⅓, ¼`), is component-dependent, and vanishes as `p → ∞`; generic `t`-dependence is `≲ 1e-3`; `o(1)` is `O(1/N)` in the last partial quotient. Evidence: V15–V21.

**C14″ (structure of Ĝ and of κ̂).** `κ̂(x̃) = lim κ` (along sequences with `p → ∞`) exists for every irrational `x̃` and every one-sided approach to a rational; `Re κ̂` is continuous with cusps at rationals of depth `≍ q'^{-2}` and Hölder exponent `≈ ½`; `Im κ̂` is odd, with **jumps** `≈ −0.05/q'` at every `p'/q'` (V23); `Ĝ` is continuous (jumps `≲ 0.002`) with `q'^{-2}` cusps. Proposed Fourier form `κ̂(x̃) = Σ_{m≥1} W_m e^{2πimx̃}/m^s + c.c.`-type with `s ∈ {1,2}` mixing Bernoulli (`B₂`, cusps) and sawtooth/Clausen (jumps) parts — the `m = kq'` terms produce the Ford weights. Open: exact cusp/jump laws in `q'` (V7 gave `≈ 0.18/q'²` for `Ĝ` on four points).

**C16′ (two-ended resonance structure, replaces C16).** `ι_{p/q}/q → Σ_{j} w(j)/(1 − λ₀^j)`-type sum over near-resonances `‖jp/q‖`. Two families of `j` are `O(1/q)`-resonant: `j ≡ mp̄ (mod q)` (`‖jp/q‖ = m/q`; gives `Ĝ(x̃)`, universal, "renormalised") and `j = k` small when `p` is bounded (`‖kp/q‖ = kp/q`; gives `Ĥ`, component-dependent because the first iterates see the component's own geometry). For `p = 1` both families coincide (explains why `p=1` is the extreme case: `G(1/q) → 1.0265` vs `≈1.017` for `[0;a₁,N]`, `a₁ → ∞`).

**C17′ (analytic mechanism, now half proven).** ~~`ρ = 1 − q²ε + κ q⁴ε² + …` with `κ` a small-denominator sum~~ → **P2**: `κ = (ι − ½)/q`. Remaining: the asymptotics `ι_{p/q} = ½ + q κ̂(x̃) + O(1)` as `q → ∞`, `p → ∞`. Route: index formula on the `W = w^q` quotient map `Φ(W) = W h(W)^q` (`ι = [W²]`-data of `Φ`), with `h` from the resonant normal form of `f_{λ₀}` whose coefficients are products of `1/(1 − λ₀^j)`; equivalently, `ι = −Σ_{other cycles of period | q} 1/(1−ρ)` [Mil] — a sum over the repelling cycles at the root, dominated by the least repelling ones. Cross-reference [BE02] for an a-priori bound `|Re ι| ≲ q`.

**C20 (`p=1` constant).** `κ₀ := lim_{q→∞} (ι_{1/q} − ½)/q = 0.0238260 − 0.052305i` (V24; approach `c/q`, `c ≈ 0.155i`) is an Écalle–Voronin / horn-map invariant of `z + z²` (parabolic implosion at `λ₀ → 1`), hence expressible via the Fatou coordinate of `z+z²`. Six digits available for identification; the naive candidates (`1/42`, `π/60`) are excluded.

**C22 (integrality).** `q·ι_{p/q}·a(p/q)² ∈ ℤ[ζ_q]` for all `q`, with denominator of `ι a²` equal to the odd part of `q`. Evidence: P3 (`q ≤ 8`); proof sketch there. Consequence: `ι` is determined by the *algebraic integer* `qιa²` of `ℚ(ζ_q)` and its Galois orbit — the object for the `q → ∞` asymptotics is a sequence of algebraic integers of degree `φ(q)`.

**C21′ (mechanism for C21, exact form).** For `λ` not a root of unity, `f_λ` has the formal linearisation `Φ_λ(z) = z + Σ φ_k(λ) z^k`, `φ_k = N_k(λ)/(λ^k − λ)` with `N_k` a finite tree sum over products of `1/(λ^{j} − λ)`, `j < k`. From `f^q = Φ^{-1}(λ^q Φ)`: `a(p/q) = [z^{q+1}] f_{λ₀}^q = lim_{λ→λ₀} (λ^q − λ^{q(q+1)}) φ_{q+1}(λ) = −q λ₀^{-1} N_{q+1}(λ₀)`. The small denominators in `N_{q+1}(λ₀)` are `1/(λ₀^{j} − 1)`, `1 ≤ j ≤ q−1`, smallest at `j ≡ ±mp̄` — the `x̃`-family; their product structure is what Yoccoz's estimate turns into the Brjuno sum (V29). The same tree sum at order `2q+1`, reduced by the index (P2 remark (i)), is `ι`; the near-resonant subsum is the candidate for `κ̂(x̃)`, the small-`j` subsum for the bounded-`p` part. This is the concrete form of next move 1.

**C21 (Yoccoz–Brjuno for the parabolic coefficient).** `log|a(p/q)| = −c·B_fin(p/q)·(1+o(1)) + O(q)`-type law with `c` universal, i.e. the leading parabolic coefficient of `e^{2πip/q}z + z²` obeys the rational analogue of Yoccoz's inequality for the Siegel radius; and `Ĝ`'s `½`-Hölder regularity (V28) is the [MMY] exponent of `log r + B`. Test: `scripts/brjuno_test.py` at `q = 127, 251` (`scripts/leading_coeff.py` first); then the exact rational version of [Yoc95] (the parabolic coefficient is the residue of the linearisation coefficient `φ_{q+1}(λ)` at `λ₀`, a small-divisor sum).

**C2′ (bulb zeta), C9′ (limit measure), C8 (Christoffel discrepancy), C11 (quadratic irrationals), C15′ (family universality):** unchanged; in C2′/C9′ read `⟨Ĝ⟩` as `⟨G⟩` (the bounded-`p` part has density zero in the Farey sequence, so the constants are unaffected); C11's rate is now seen for `κ` (V20). C15′'s question "is the cusp law family-universal?" is sharpened by V19: for `z³+c` the `x̃`-part should be a *different* universal function `Ĝ₃`, with its own bounded-`p` part.

---

## 6. Relationship to Ford circles (current statement)

| level | Ford circles | bulbs |
|---|---|---|
| combinatorics | Farey adjacency, mediants | identical |
| leading order | radius `1/(2q²)`, exact | `2|c_H'|q⁻²`, asymptotic [GM84] |
| first correction | none (PSL(2,ℤ)-exact) | `ρ = 1 − q²ε + q³(ι−½)ε²`, `ι ∈ ℚ(ζ_q)` a Galois orbit; `κ̂(p̄/q)` singular on ℚ with weights decaying in `q'` |
| symmetry | PSL(2,ℤ) on horoballs | Galois `ζ_q ↦ ζ_q^{p'p̄}` on `ι`; `p ↦ p̄` as the dominant variable, `p ↔ p̄` the two resonance ends |

---

## 7. Next moves (ordered)

1. **C17′ → asymptotics of `ι`.** The cycle-index sum is collective (V22), so the route is the normal form: `ι = [W²]`-data of `Φ(W) = W h(W)^q` with `h` from the resonant normal form of `f_{λ₀}`, whose coefficients are products of `1/(1 − λ₀^j)`; extract the `j ≡ mp̄` (near-resonant) and small-`j` (bounded-`p`) families and compare with `κ̂(x̃)` and `Ĥ(p)` from V16/V18. First test: does the *first* resonant coefficient `b₁` alone (`|b₁|`, `arg b₁` vs `x̃`) already carry the `x̃`-structure? Cheap with `scripts/exact_index.py` machinery.
2. **Cusp profile → tip function.** Settled: V-shaped along each family (V30), jumps `J ≈ −0.05/q'` in `Im κ̂` (V23). Open: the *tip map* `a ↦ G_tip(p'/q'; a)` (the deep-digit dependence, `≈ 0.001–0.004`) and the family slopes `s(p'/q')` (`1.65, 1.29, 0.58` at `⅓, ½, 0`): are they `q'^{-1}`-, `q'^{-2}`- or Brjuno-weighted? Cheap: repeat `canonical_profile.py` at `¼, ⅕, ⅖, 2/7` and with `a ∈ {2, 3, 5, 16, 64}`.
3. **Identify `κ₀`** (C20): Fatou-coordinate / horn-map computation for `z + z²`, matched against the six digits of V24.
4. **`z³+c`**: repeat V15–V19 (`MAIN3` is already wired: `G(MAIN3, p, q)`, `taylor(p, q, MAIN3)`).
5. C2′ numerics, then C8, C11 as before.

## 8. Instruments (`bulbford/`, tests in `tests/`, 17 passing)

- `cf.py` — `modinv`, `xstar`, `cf`, `from_cf`, `convergent_denominators`.
- `dynamics.py` — `Family` records (`MAIN2`, `DISK2`, `MAIN3`), `orbit` with second derivatives, `Cycle` tracking with analytic `dρ/dc` and a period-collapse guard, `bulb(fam, p, q)` → `G_ant`, `G_cen`; `rho_on_path`.
- `index.py` — `index(p, q)`: `ι_{p/q}` as an `acb` ball (python-flint), auto precision; `kappa(p, q)`.
- `taylor.py` — `taylor(p, q, fam, r, N)`: `r_k`, `solve(target, u0)`; `kappa_fft` (N=64) for `q > 300`.
- `scripts/` — `exact_index.py` (ℚ(ζ_q)), `taylor_sweep.py`, `index_sweep.py`, `kappa_sweep.py` (dense, `q = 1009, 2003`), `sequences.py`, `prefix_test.py`, `a1_scan*.py`, `limit_grid.py`, `generic_limits.py`, `p1_limit.py`, `cycle_decomposition.py`, `cusp_fit.py`, `main3_limits.py`, `analyze.py`, `analyze_dense.py`. Data in `data/` (JSON). Legacy instruments in `legacy/`.
- Run: `PYTHONPATH=. pytest -q`; `PYTHONPATH=. python3 scripts/analyze.py 59 127 251`.

## 9. Error / caveat catalog

- E-a. Diameter proxy (ρ=−1 point). Unchanged; `|u_a|/2` is the cleaner proxy (V14).
- E-b. V7 rests on four cusp depths; the jump/cusp laws are open (next move 2).
- E-c. `q > 2100` untested; period-collapse guard bisects continuation steps (no failures logged).
- E-d. ~~`Ĝ(½)` extrapolated~~ → direct limit `1.11142` along `p=2`; other approaches to `½` differ at the V17 level.
- E-e. "Converged" for sequences means stable to `1e-5` across a doubling of `N`; the `[0;…,N,a]` families converge only like `c/N` and are quoted with the fit.
- E-f. The `ι`-route needs `dps ≈ 1.5q` and `O(q³)` ball operations; used to `q = 251` (385 s). Beyond that the FFT route is the instrument (validated to `1e-12` by V12).
