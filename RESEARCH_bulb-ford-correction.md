# RESEARCH BRANCH — Bulb geometry of M vs Ford circles: the arithmetic correction Ĝ

Branch id: `bulb-ford-correction`  ·  opened 2026-09-24  ·  revised 2026-09-25, 2026-09-27  ·  status: **active, theorems P2 and P6, numerics-validated to q ≈ 2000**
Limit-free companion: [`RESEARCH_bulb-ford-correction_finite.md`](RESEARCH_bulb-ford-correction_finite.md), the same claims with no limit asserted (enforced by `tools/audit_limits.py`). The 2026-09-27 revision folds its corrections back into this file: the P2 proof, P6, V14, F4 and F5.
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

Theorem tags: **[GM84]** Guckenheimer–McGehee `diam B_{p/q} = O(q⁻²)`; **[DH]** Douady–Hubbard multiplier map; **[BS94]** Bullett–Sentenac; **[Yoc95]** Yoccoz; **[MMY]** Marmi–Moussa–Yoccoz; **[Weil]** Kloosterman; **[FL]** Franel–Landau; **[Mil]** Milnor, *Dynamics in One Complex Variable*, §12 (holomorphic index, `ι = b/a²` for `z + az² + bz³`, index sum `Σ 1/(1−ρ) = 0` over finite fixed points of a polynomial); **[Gol92]** Goldberg, rotation sets of `θ ↦ 2θ`; **[Mil00]** Milnor, *Periodic orbits, external rays and the Mandelbrot set*; **[BE02]** Buff–Epstein, parabolic Pommerenke–Levin–Yoccoz inequality (bounds `Re resit` for quadratic parabolic points — *statement to be re-checked before use*).

---

## 2. PROVEN (this branch)

- **P0.** First-order local model: `ρ = 1 − q(λ^q − 1) + O((λ^q−1)²) = 1 − q²ε + O(ε²)`; leading law `d = 2|c_H'(λ₀)| q⁻²` to first order. (Elementary.)
- **P1.** `q = 2` exact: period-2 disk diameter `1/2 = 2·sin(π/2)/4`. Also `ρ(λ) = −λ² + 2λ + 4` for the 2-cycle of `λz+z²` (tested).
- **P2 (index formula for the second-order term).** For any unicritical family and any hyperbolic component `H`, with `λ = λ₀e^ε`:
  `ρ(ε) = 1 − q²ε + q³(ι − ½) ε² + O(ε³)`, `ι = ι(f_{λ₀}^{q·per(H)}, z₀)` the holomorphic index of the parabolic cycle point.
  ~~*Proof.* `ρ(λ)` is analytic at `λ₀` (symmetric function of the `q` cycle points, which are the roots of a dynatomic factor inside a fixed disc; residue-theorem continuity). The `q+1` fixed points of `F = f^q` merging at `z₀` — the `H`-cycle point (multiplier `μ = λ^q = e^{qε}`) and the `q` satellite points (multiplier `ρ`) — have total index converging to the index of the merged parabolic point [Mil, 12.9]: `1/(1−μ) + q/(1−ρ) → ι`. With `ρ = 1 − q²ε + Cε² + O(ε³)`: `1/(1−μ) = −1/(qε) + ½ + O(ε)`, `q/(1−ρ) = 1/(qε) + C/q³ + O(ε)`, hence `½ + C/q³ = ι`. ∎~~ (superseded 2026-09-27: correct, but it assumes the linear coefficient and passes through a limit.)
  *Proof (revised; one Laurent identity, which also proves P0's linear term).* Fix a circle `γ` around `z₀` enclosing no other fixed point of `F₀`, and put `I(ε) = (1/2πi)∮_γ dz/(z − F_ε(z))`. Then `I` is holomorphic in `ε` with `I(0) = ι`, and `ρ` is analytic at `ε = 0` (symmetric function of the cycle points). The argument principle gives exactly `q + 1` fixed points inside `γ`, which by [DH] are the `H`-point (multiplier `μ = e^{qε}`) and the `q` satellite points (multiplier `ρ`). All are simple for small `ε ≠ 0`, so the residue theorem gives the identity of Laurent series `1/(1−μ) + q/(1−ρ) = I(ε)`. Here `1/(1−μ) = −1/(qε) + ½ − qε/12 + …`. Write `ρ = 1 + ρ_k ε^k + …` with `ρ_k ≠ 0`. The right side has no pole, and a pole of order `k ≥ 2` in `q/(1−ρ)` could not cancel, so `k = 1`. Cancelling the simple pole gives `ρ₁ = −q²`, and the constant term gives `½ + qρ₂/ρ₁² = ι`, that is `ρ₂ = q³(ι − ½)`. ∎
  Checks: `q=1`: `ρ = 2 − e^ε = 1 − ε − ε²/2`, `ι(z+z²) = 0` ✓. `q=2`: `ρ = 1 − 4ε − 3ε²`, `ι(z − 2z³ + z⁴) = 1/8` ✓. Numerically for all `p`, `q ≤ 251` (V12).
  *Remarks.* (i) The `ε`-derivatives of the normal-form coefficients cancel at second order (checked by direct normal-form expansion: `ρ = 1 − q²ε + q²(b₂/b₁² − ½)ε² + …` for `f ~ λw(1 + b₁w^q + b₂w^{2q})`, and `b₂/b₁² = qι − (q−1)/2`). (ii) `ι_{p/q} ∈ ℚ(ζ_q)` and `ι_{p'/q} = σ(ι_{p/q})` for the Galois automorphism `ζ ↦ ζ^{p'p̄}` (the residue is a rational function of `λ₀` with rational coefficients). (iii) Third and higher orders depend on the unfolding, not on the germ alone.

- **P3 (certified satellite centres, finite).** For every `p/q` with `2 ≤ q ≤ 16`, `gcd(p,q) = 1` (79 cases), a dyadic box `β_{p/q}` of half-width `2⁻⁶⁴` satisfies the joint witness of the finite-Mandelbrot certificate calculus for type `(0,q)`: `K_{Q_q}(β) ⊂ int β` (a unique, simple root of the Gleason polynomial in `β`) and `0 ∉ Q_j(β) − Q_i(β)` for every forbidden pair (exact critical period `q`), on the same box; the boxes of one `q` are pairwise disjoint. Checked in rational interval arithmetic, replayed from stored endpoints (`data/center_certificates.json`, `tests/test_certify.py`). *Imports for reading `β_{p/q}` as the centre of `B_{p/q}`*: [DH], and the continuation seed (`SatelliteLabel`, VALIDATED only).
- **P4 (wake combinatorics, finite).** For `q ≤ 16`: the doubling cycle with mechanical word `b_k = [kp mod q ≥ q−p]` acts as the rotation by `p/q`; its shortest arc `(θ₋, θ₊)` has width `1/(2^q − 1)`; `(q−p)/q` gives the mirror arc; the arcs of `F_n` (`n ≤ 12`) are disjoint and ordered like the fractions, and for Farey neighbours (`n ≤ 10`) the mediant's arc lies in their gap (`tests/test_wake.py`). *Imports for reading the arc as the `p/q`-wake*: [Gol92] (uniqueness of the rotation cycle), [DH]/[Mil00] (landing at the root). This is the "combinatorics: identical" row of §6 made executable.
- **P6 (projective index sum; added 2026-09-27).** The homogeneous lift `[Z:W] ↦ [W^d F(Z/W) : W^d]` of a polynomial `F` of degree `d ≥ 2` fixes `∞ = [1:0]` with multiplier `0`, so `ι_∞ = 1`. The index formula `Σ_{ℙ¹} ι = 1` [Mil, §12] then gives `Σ_{z∈ℂ} ι_z = 0`. For `F = f_{λ₀}^q`, every fixed point other than `z₀` is repelling [DH], hence simple, so `ι_{p/q} = −Σ_{F(z)=z, z≠z₀} 1/(1 − F'(z))`, exactly, over finitely many points. Evaluated for `q ≤ 20` in V22.
- **P5 (exact `ι`, `q ≤ 8`).** The `finite-math-kernels` cyclotomic reference (`[w^q] 1/P` for `ζ_q^p w + w²`, exact in `ℚ(ζ_q)`, every unit `p`) equals `scripts/exact_index.py` coordinate-by-coordinate, lies in the Arb ball of `index.py`, and reproduces the pinned `ℚ(i)` fixtures (`tests/test_cyclotomic_replay.py`); `ι_{p/q} = σ_p(ι_{1/q})` holds exactly (kernel C2 tests). This makes remark (ii) of P2 a checked identity for `q ≤ 8`.

---

## 3. VALIDATED (numerical; doubles unless stated; Newton to 1e-14 relative)

Legacy items V1–V11 stand (reproduced 2026-09-25 by the new instruments: V5/V9/V10 table at `q=59` identical to 4 digits). Corrections: V6's "one-sided limits agree" is superseded by V16; V4's "residual dependence on `a_{n−1}` ≈ 0.005" is explained by V17–V18.

- **V12. P2 verified.** `[u²]R_q` (Cauchy/FFT on `|u| = 2.5`, 256 points) vs `(ι−½)/q` (ball arithmetic, radius `< 1e-12`): max deviation over all `p`: `9e-15 (q=59)`, `1e-13 (127)`, `1.3e-12 (251)`. `r₀ = 1`, `r₁ = −1` to `1e-11`.
- **V13. Exact values** (`q ≤ 12`, `data/exact_index.txt`): `ι_{1/3} = (92 − 16ζ₃)/441`, `ι_{1/4} = (1447 − 365i)/4624`, `ι_{1/5} = (1068·31 − 9153ζ − 7113ζ² − 312ζ³)/93775`; denominators are squares of algebraic norms (`441 = 21²`, `4624 = 68²`, `3025 = 55²`), as expected for `b/a²`. `Im ι(1/q)` decreases by `0.064 ± 0.001` per unit `q` for `q = 3…12`: linear growth, i.e. `κ = O(1)` already visible exactly.
- **V14. The limit shape lives in the `ε`-plane.** `|u_a|/2` (root of the truncated Taylor series of `R_q` at `u=0`, radius of the sampling circle 2.5) reproduces `G_ant` (direct antipode continuation) to `5.7e-4 (q=59)`, `1.3e-4 (127)`, `5.5e-5 (251)` — ~~i.e. `O(q⁻²)`, the nonlinearity of `c_H(λ)`~~ (withdrawn 2026-09-27: `q²Δ = 1.98, 2.10, 3.47` is not constant, so these three rows do not establish the rate). `|r_k|` decays like `q^{-k}`-ish (tail `< 1e-16` at `k=16`); `R_q` is analytic on `|u| ≲ cq`.
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
  Dense `q = 1009` (all 504 `p ≤ q/2`, `data/kappa_q1009.json`): `var G = 1.64e-3`; a 100-bin function of `x*` alone leaves residual variance `2.65e-6` (`R² = 0.9984`, residual sd `1.6e-3`); adding a 50-bin function of `t` (backfitting) leaves `1.87e-6` (`R² = 0.9989`, `sd H = 1.0e-3`). So generically the `t`-part is `≈ 1e-3`. The large residuals (local: `G` minus its 8 nearest `x*`-neighbours) are the bounded-`p` bulbs: `p = 1,2,3,4`: `+0.0080, +0.0044, +0.0025, +0.0005`; and the bulbs with `p·p̄ = q−1` (`x̃ = −1/p` exactly, `6 ≤ p ≤ 36`): `−0.0036 … −0.0053` (they sit at the bottom of the `1/p` cusps). Small-denominator `t` classes (`t ≈ ⅓, ¼, ⅖`, `n = 4` each) average `−0.002 … −0.003` (2–3σ). Generic `p ≥ 20` off these sets: `|e| < 1e-3`.
- **V21. `Im κ` jumps at rationals `p'/q'` (q = 1009, means over `0.002 < |x̃ − p'/q'| < 0.009` on each side, `J := Im κ(p'/q'⁺) − Im κ(p'/q'⁻)`):**
  `½: −0.034 (from the two conjugate sides)`, `⅓: −0.0186`, `¼: −0.0105`, `⅕: −0.0045`, `⅖: −0.0053`, `⅙: −0.0029`, `⅐: −0.0004`, `2/7: −0.0007`, `3/7: 0.0000`.
  `J·q'² = −0.137, −0.167, −0.167, −0.112, −0.132, −0.106` for `q' = 2,3,4,5,5,6`: Ford-type `q'^{-2}` weight, equal for `⅕` and `⅖`, then falling faster at `q' = 7` (window-limited at this `q`). Same sign everywhere (the jump is *down* in `Im κ` for `x̃ > 0`; odd in `x̃`). Together with V16(b): `Re κ` has the cusps, `Im κ` the jumps.
- **V19. Universality across components is exact for the `x̃`-part and fails only for bounded `p`.** `disk2/main2` ratio of `G` on `[0;a₁,16,3]`: `0.9924` at `p=3` (`a₁=1`), `1.0000 ± 0.0002` for `a₁ = 2…100`; on `[0;a₁,64]`: `0.9881` at `p=1`, `0.9997–1.0008` for `a₁ ≥ 2`. Refines V9 (`0.9985 ± 0.0032` at `q=59`, whose scatter was the bounded-`p` bulbs).
- **V20. Quadratic irrational:** `x̃ → 2−φ = [0;2,1,1,…]`: `κ = 0.05716, 0.06216, 0.06390` at `q = 34, 89, 610`; geometric convergence (ratio `≈ 1/3` per two Fibonacci steps), `κ_∞ ≈ 0.0648 − 0.0058i`, `G_∞ ≈ 1.158` (the global maximum region).
- **V22. P6 evaluated (`q ≤ 20`, all `p`; `bulbford/cycles.py`, `data/cycle_index_sum.json`).** `ι_{p/q} = −Σ 1/(1 − F'(z))` over the `2^q − q − 1` other fixed points of `f^q`, one per ray of period dividing `q` outside the `p/q` rotation cycle, agrees with the Arb index to `≤ 1.8e-14` on all 64 fractions with `p ≤ q/2`. **Flank cycle:** let `(α_w, α_{w+1})` be the characteristic arc of the rotation-cycle angles `α_0 < … < α_{q−1}`. The doubling orbit of `α_{w−1} + 1/M` passes through `α_{w+2} − 1/M` (exact in integers, every `p`, `3 ≤ q ≤ 20`), and this period-`q` cycle (not a rotation cycle for `q ≥ 4`) carries the largest term in all 63 fractions with `q ≥ 3`. Its share is `0.04–0.11` at `q ≥ 17`, and `|F'| ∈ [14.9, 51.3]`. For `p ≥ 2` at fixed `q`, the flank `|F'|` is strictly decreasing in `x*` (30/30 comparisons, `q = 5…20`), so the dominant term is organised by the same variable as `Ĝ`. `p = 1` is off the order (bounded-`p`, cf. V16). For `p = 1`, `|F'|` = `9.00, 14.93, …, 48.31` for `q = 2…20`, increasing and flattening. The sum cancels heavily (`Σ|term|/|ι|` up to `7.4`; the top ten terms carry `0.12–0.20` for `p ≥ 2` at `q ≥ 19`). Rotation cycles carry `≤ 0.084` (`q ≥ 11`), and `β` carries `≤ 0.047` (`q ≥ 6`).
- **V23. Flank cycle against `κ` and `G` (`q ≤ 20`; `scripts/flank_vs_kappa.py`, `data/flank_vs_kappa.json`).** Over the 78 within-`q` pairs of fractions with `p ≥ 2` (`q = 5…20`), the flank `|F'|` is anti-ordered with `Re κ` on 67 pairs and with `G_ant` on 68. Every exception involves `p = 2`, the bulb with `x*` nearest `½`, where `Re κ` and `G` have passed their peak (V15: maximum near `x* ≈ 0.40`) while the flank `|F'|` keeps falling (strictly decreasing in `x*`, V22(d)). Example, `q = 19`: `x* = 0.105 … 0.421` gives flank `51.29 → 34.91`, `Re κ` `0.0309 → 0.0490` and `G` `1.0560 → 1.1188`. At `x* = 0.474` the flank is `33.92` but `Re κ = 0.0443` and `G = 1.1031`. **Reading:** the dominant term of the P6 sum is a monotone function of `x*`, and `Re κ`, being unimodal, is not a function of the flank term alone. The flank term orders the rising branch of `Ĝ` but not its fall toward `½`.

Precision caveats: doubles for `R_q`; the FFT route is accurate to `~1e-12` in `r₂` (verified against balls, V12) and is the instrument for `q > 300`. `G` from `|u_a|/2` carries the `O(q⁻²)` proxy difference of V14 (irrelevant for limits). Prime `q` used in sweeps; sequences use whatever `q` the CF gives.

---

## 4. FALSIFIED (this branch)

- F1 (was C12). `c₁(t)/q` correction with smooth `c₁`. **False** (V2–V3).
- F2 (was C13). "Nothing arithmetic transfers." **False** (V3–V6, V15).
- F3 (was C15). Same `Ĝ` for every degree. **False** (V10).
- **F4 (was C1″/C14′ as stated).** "`d = 2|φ'|q⁻² Ĝ(x*)(1+o(1))` with a single continuous function `Ĝ` on `(0,½]`." **False at the 0.5–1.8 % level**: the `q→∞` limit along `x̃ → x` depends on whether `p` stays bounded (V16–V17: `1.0265` vs `≈1.008` at `0⁺`, `1.1284` vs `1.1205` at `⅓⁻`), and one-sided generic limits at `⅓` differ by `≈0.002` (jump on top of the cusp). The leading `x*`-description survives as the dominant term (range `0.15` vs generic residual `1.6e-3`, R² = 0.9984 at `q = 1009`). *Finite witness (2026-09-27):* bulb pairs with `|Δx̃| ≈ 1e-7` differ in `G` by up to `0.016` (`1/256`: `1.02628` vs `16129/16385`: `1.00998`). Any `g` with `G = g(x̃)` on the 678 computed bulbs therefore has Lipschitz constant `≥ 6.8·10⁴` (table in the finite register, F4).
- **F5 (was the contraction form of C16).** "Deeper ancestors enter with weight `(q_{n−1}/q)²`." **False**: the digit `a₁` of `p/q` (deepest ancestor of `x̃`) changes the limit by `0.005–0.012` even as the intermediate digit `N → ∞` (V17), which would be suppressed by `(1/N)²` under C16. *Finite form (2026-09-27):* the model's constant must satisfy `C ≥ gap·N² = 1.8, 20.8, 108` at `N = 32, 64, 128` (V17 grid), which grows with `N`.
- **F6 (the Farey-neighbour expectation of next move 1).** "The cycles of the Farey neighbours carry the largest share of `ι`." **False as stated**: for prime `q`, neighbours with `1 < q' < q` have `q' ∤ q` and no term in the sum, the `q' = 1` neighbour `β` carries `≤ 0.047` for `q ≥ 6`, and the largest term is the flank cycle in every case (V22).
- ~~V6's claim of equal one-sided limits at ⅓~~ (superseded; the `0.001` agreement there was at one fixed prefix and inside its error).

---

## 5. CONJECTURED (open register)

**C1‴ (size law, corrected form).** For every hyperbolic component `H` of the quadratic family,
`d_H(p/q) = 2|c_H'(λ₀)| q⁻² · G_H(p/q) · (1 + O(q⁻²))`, with
`G_H(p/q) = Ĝ(x̃) + Ĥ_H(p; x̃) + o(1)`, where `Ĝ` is universal (component-independent), depends on `x̃ = p̄/q` only, range `≈[1.01, 1.16]`, and is the limit along any sequence with `p → ∞`; `Ĥ_H` is supported on bounded `p` (`+0.018, +0.011, +0.008, +0.0005` for `p = 1,2,3,4` at `x̃ → 0, ½, ⅓, ¼`), is component-dependent, and vanishes as `p → ∞`; generic `t`-dependence is `≲ 1e-3`; `o(1)` is `O(1/N)` in the last partial quotient. Evidence: V15–V21.

**C14″ (structure of Ĝ and of κ̂).** `κ̂(x̃) = lim κ` (along sequences with `p → ∞`) exists for every irrational `x̃` and every one-sided approach to a rational; `Re κ̂` is continuous on irrationals with V-cusps at rationals; `Im κ̂` is odd, with **jumps** at rationals `p'/q'`; both singular amplitudes decay in `q'`. Proposed Fourier form `κ̂(x̃) = Σ_{m≥1} W_m e^{2πimx̃}/m^s + c.c.`-type with `s ∈ {1,2}` mixing Bernoulli (`B₂`, cusps) and sawtooth/Clausen (jumps) parts — the `m = kq'` terms produce the Ford weights. Open: exact cusp/jump laws in `q'` (V7 gave `≈ 0.18/q'²` for `Ĝ` on four points).

**C16′ (two-ended resonance structure, replaces C16).** `ι_{p/q}/q → Σ_{j} w(j)/(1 − λ₀^j)`-type sum over near-resonances `‖jp/q‖`. Two families of `j` are `O(1/q)`-resonant: `j ≡ mp̄ (mod q)` (`‖jp/q‖ = m/q`; gives `Ĝ(x̃)`, universal, "renormalised") and `j = k` small when `p` is bounded (`‖kp/q‖ = kp/q`; gives `Ĥ`, component-dependent because the first iterates see the component's own geometry). For `p = 1` both families coincide (explains why `p=1` is the extreme case: `G(1/q) → 1.0265` vs `≈1.017` for `[0;a₁,N]`, `a₁ → ∞`).

**C17′ (analytic mechanism, now half proven).** ~~`ρ = 1 − q²ε + κ q⁴ε² + …` with `κ` a small-denominator sum~~ → **P2**: `κ = (ι − ½)/q`. Remaining: the asymptotics `ι_{p/q} = ½ + q κ̂(x̃) + O(1)` as `q → ∞`, `p → ∞`. Route: index formula on the `W = w^q` quotient map `Φ(W) = W h(W)^q` (`ι = [W²]`-data of `Φ`), with `h` from the resonant normal form of `f_{λ₀}` whose coefficients are products of `1/(1 − λ₀^j)`; equivalently, `ι = −Σ_{other cycles of period | q} 1/(1−ρ)` [Mil] — a sum over the repelling cycles at the root, dominated by the least repelling ones. Cross-reference [BE02] for an a-priori bound `|Re ι| ≲ q`.

**C20 (`p=1` constant).** `κ₀ := lim_{q→∞} (ι_{1/q} − ½)/q = 0.0238 − 0.0527i` (Im to `±0.001`) is an Écalle–Voronin / horn-map invariant of `z + z²` (parabolic implosion at `λ₀ → 1`: the `1/q`-cycle's multiplier is governed by the Lavaurs phase), hence expressible via the Fatou coordinate of `z+z²`. Untested.

**C2′ (bulb zeta), C9′ (limit measure), C8 (Christoffel discrepancy), C11 (quadratic irrationals), C15′ (family universality):** unchanged; in C2′/C9′ read `⟨Ĝ⟩` as `⟨G⟩` (the bounded-`p` part has density zero in the Farey sequence, so the constants are unaffected); C11's rate is now seen for `κ` (V20). C15′'s question "is the cusp law family-universal?" is sharpened by V19: for `z³+c` the `x̃`-part should be a *different* universal function `Ĝ₃`, with its own bounded-`p` part.

---

## 6. Relationship to Ford circles (current statement)

| level | Ford circles | bulbs |
|---|---|---|
| combinatorics | Farey adjacency, mediants | identical (P4: exact wake arcs, Farey order, mediant in the gap) |
| leading order | radius `1/(2q²)`, exact | `2|c_H'|q⁻²`, asymptotic [GM84]; on every computed bulb `G ∈ [1.0144, 1.1571]` (`q = 1009`) |
| first correction | none (PSL(2,ℤ)-exact) | `ρ = 1 − q²ε + q³(ι−½)ε²`, `ι ∈ ℚ(ζ_q)` a Galois orbit; `κ̂(p̄/q)` singular on ℚ with weights decaying in `q'` |
| symmetry | PSL(2,ℤ) on horoballs | Galois `ζ_q ↦ ζ_q^{p'p̄}` on `ι`; `p ↦ p̄` as the dominant variable, `p ↔ p̄` the two resonance ends |

---

## 7. Next moves (ordered)

1. **C17′ → asymptotics of `ι`.** ~~Farey-neighbour cycles dominate~~ (F6). Done for `q ≤ 20` (V22, V23): the flank cycle dominates, its `|F'|` is monotone in `x*`, and it orders the rising branch of `Re κ` and `G` but not their fall toward `½`. Next: prove the flank-orbit identity in `ℚ/ℤ` for all `q`, identify the second-rank terms (one `α`-neighbour each, gap pattern `1, 2, 4, 8`), and push the ray method to `q = 24` (16M rays).
2. **Cusp/jump laws.** From `data/kappa_q1009.json`, fit `Re κ̂` cusps and `Im κ̂` jumps at `p'/q'` for `q' ≤ 7` against `q'^{-1}`, `q'^{-2}`; confirm with `q = 2003`.
3. **Fatou-coordinate computation of `κ₀`** (C20), cheap via the standard `z + z²` Fatou coordinate series.
4. **`z³+c`**: repeat V15–V19 (`MAIN3` is already wired: `G(MAIN3, p, q)`, `taylor(p, q, MAIN3)`).
5. C2′ numerics, then C8, C11 as before.
6. **Certified antipodes.** Extend P3 to `c_ant` (`ρ = −1`) by a two-variable Krawczyk witness for `(f_c^q(z) − z, (f_c^q)'(z) + 1)` on a box in `(z, c)`: with [DH] this turns `G_ant` for small `q` into a certified enclosure. Then the exact root `c_root = λ₀/2 − λ₀²/4 ∈ ℚ(ζ_q)` via a certified `ζ_q` box (Krawczyk on `Φ_q`, no trigonometry).

## 8. Instruments (`bulbford/`, tests in `tests/`, 381 passing)

- `cf.py` — `modinv`, `xstar`, `cf`, `from_cf`, `convergent_denominators`.
- `dynamics.py` — `Family` records (`MAIN2`, `DISK2`, `MAIN3`), `orbit` with second derivatives, `Cycle` tracking with analytic `dρ/dc` and a period-collapse guard, `bulb(fam, p, q)` → `G_ant`, `G_cen`; `rho_on_path`.
- `index.py` — `index(p, q)`: `ι_{p/q}` as an `acb` ball (python-flint), auto precision; `kappa(p, q)`.
- `certify.py` — rational dyadic interval boxes, `Q_n`/`Q_n'` jets on boxes, collision sets `I/F_{ℓ,k}(H)`, Krawczyk, joint certificate with verdict `ACCEPTED`/`INCONCLUSIVE`, record replay; `scripts/certify_centers.py` (`--check`).
- `wake.py` — rotation cycles of doubling, characteristic arcs, Farey fractions.
- `taylor.py` — `taylor(p, q, fam, r, N)`: `r_k`, `solve(target, u0)`; `kappa_fft` (N=64) for `q > 300`.
- `cycles.py` — P6: all fixed points of `f^q` by vectorized ray pull-back, `cycle_terms`, `index_by_cycles`, `flank_angles`; `scripts/cycle_index_sweep.py`, `scripts/flank_vs_kappa.py`.
- `tools/audit_limits.py` — the no-limits audit of the finite register (CI).
- `scripts/` — `exact_index.py` (ℚ(ζ_q)), `taylor_sweep.py`, `index_sweep.py`, `kappa_sweep.py` (dense, `q=1009`), `sequences.py`, `prefix_test.py`, `a1_scan*.py`, `limit_grid.py`, `analyze.py`, `analyze_dense.py`. Data in `data/` (JSON). Legacy instruments in `legacy/`.
- Run: `PYTHONPATH=. pytest -q`; `PYTHONPATH=. python3 scripts/analyze.py 59 127 251`.

## 9. Error / caveat catalog

- E-a. Diameter proxy (ρ=−1 point). Unchanged; `|u_a|/2` is the cleaner proxy (V14).
- E-b. V7 rests on four cusp depths; the jump/cusp laws are open (next move 2).
- E-c. `q > 2100` untested; period-collapse guard bisects continuation steps (no failures logged).
- E-d. ~~`Ĝ(½)` extrapolated~~ → direct limit `1.11142` along `p=2`; other approaches to `½` differ at the V17 level.
- E-e. "Converged" for sequences means stable to `1e-5` across a doubling of `N`; the `[0;…,N,a]` families converge only like `c/N` and are quoted with the fit.
- E-f. The `ι`-route needs `dps ≈ 1.5q` and `O(q³)` ball operations; used to `q = 251` (385 s). Beyond that the FFT route is the instrument (validated to `1e-12` by V12).
