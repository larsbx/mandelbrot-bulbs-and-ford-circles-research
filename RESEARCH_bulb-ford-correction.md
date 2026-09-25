# RESEARCH BRANCH — Bulb geometry of M vs Ford circles: the arithmetic correction Ĝ

Branch id: `bulb-ford-correction`  ·  opened 2026-09-24  ·  revised 2026-09-25  ·  status: **active, one theorem (P2), numerics-validated to q ≈ 2000**
Discipline: PROVEN / VALIDATED (numerical, stated precision) / CONJECTURED / FALSIFIED are absolute. Analytic imports carry theorem tags. Superseded items are struck, never deleted.

---

## 0. One-paragraph summary

The satellite bulbs of every hyperbolic component of the quadratic family obey
`d(p/q) = 2|φ_H'(e^{2πip/q})| · q⁻² · G(p/q) · (1+O(q⁻²))`, with `G` an `O(1)` correction whose leading dependence is on the modular inverse `x̃ = p̄/q ∈ (−½,½]`, `p·p̄ ≡ 1 (mod q)`: range 0.12 in `|x̃|`, a downward cusp at every rational. **New (P2):** the whole second-order term is controlled by one classical invariant — the holomorphic index `ι_{p/q}` of the parabolic point of `e^{2πip/q}z + z²` — through `ρ = 1 − q²ε + q³(ι − ½)ε² + O(ε³)`; `κ := (ι − ½)/q` is `O(1)`, and `G` is a function of `κ` and the next Taylor coefficient to `10⁻³`. **Correction to the opening statement:** `G` is *not* a function of `x*` alone. A second, weaker (`≤ 0.013`) dependence sits at the *other* end of the continued fraction (on `t = p/q`, through bounded `p` / small-denominator `t`), has the same arithmetic (peaks at rationals) structure, and is the *only* part that is not universal across hyperbolic components. The Ford-circle analogy is exact at the combinatorial level, asymptotic at leading order, and at first correction the modular group enters through `p ↦ p⁻¹ (mod q)` in the arithmetic of a Galois orbit in `ℚ(ζ_q)`.

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
  *Remarks.* (i) The `ε`-derivatives of the normal-form coefficients cancel at second order (checked by direct normal-form expansion: `ρ = 1 − q²ε + q²(b₂/b₁² − ½)ε² + …` for `f ~ λw(1 + b₁w^q + b₂w^{2q})`, and `b₂/b₁² = qι − (q−1)/2`). (ii) `ι_{p/q} ∈ ℚ(ζ_q)` and `ι_{p'/q} = σ(ι_{p/q})` for the Galois automorphism `ζ ↦ ζ^{p'p̄}` (the residue is a rational function of `λ₀` with rational coefficients). (iii) Third and higher orders depend on the unfolding, not on the germ alone.

---

## 3. VALIDATED (numerical; doubles unless stated; Newton to 1e-14 relative)

Legacy items V1–V11 stand (reproduced 2026-09-25 by the new instruments: V5/V9/V10 table at `q=59` identical to 4 digits). Corrections: V6's "one-sided limits agree" is superseded by V16; V4's "residual dependence on `a_{n−1}` ≈ 0.005" is explained by V17–V18.

- **V12. P2 verified.** `[u²]R_q` (Cauchy/FFT on `|u| = 2.5`, 256 points) vs `(ι−½)/q` (ball arithmetic, radius `< 1e-12`): max deviation over all `p`: `9e-15 (q=59)`, `1e-13 (127)`, `1.3e-12 (251)`. `r₀ = 1`, `r₁ = −1` to `1e-11`.
- **V13. Exact values** (`q ≤ 12`, `data/exact_index.txt`): `ι_{1/3} = (92 − 16ζ₃)/441`, `ι_{1/4} = (1447 − 365i)/4624`, `ι_{1/5} = (1068·31 − 9153ζ − 7113ζ² − 312ζ³)/93775`; denominators are squares of algebraic norms (`441 = 21²`, `4624 = 68²`, `3025 = 55²`), as expected for `b/a²`. `Im ι(1/q)` decreases by `0.064 ± 0.001` per unit `q` for `q = 3…12`: linear growth, i.e. `κ = O(1)` already visible exactly.
- **V14. The limit shape lives in the `ε`-plane.** `|u_a|/2` (root of the truncated Taylor series of `R_q` at `u=0`, radius of the sampling circle 2.5) reproduces `G_ant` (direct antipode continuation) to `5.7e-4 (q=59)`, `1.3e-4 (127)`, `5.5e-5 (251)` — i.e. `O(q⁻²)`, the nonlinearity of `c_H(λ)`. `|r_k|` decays like `q^{-k}`-ish (tail `< 1e-16` at `k=16`); `R_q` is analytic on `|u| ≲ cq`.
- **V15. `κ` is `O(1)`, and `G` is a function of `(κ, r₃)`.** All `p`, `q ∈ {59,127,251,1009}`: `Re κ ∈ [0.015, 0.064]`, `|Im κ| ≤ 0.055`, `|r₃| ≤ 2.3e-3`, `|r₄| ≤ 1e-4`. `Re κ` is a smooth unimodal function of `x*` with the profile of `Ĝ` (`0.024` at `x*→0`, max `≈0.062` near `0.40`, `0.050` at `½`). `Im κ` is **odd in the signed inverse** `x̃` (`κ(−x̃) = conj κ(x̃)`, exact by conjugation), `|Im κ| ≈ 0.053` at `x*→0`, crossing `0` near `x* ≈ 0.32–0.39`, `0.024` at `½`. Newton on `1 − u + κu² + r₃u³ = −1` from `u=2` gives `G` to `≤ 3e-3` (peak `1.147` ✓, `x*→0`: `1.029` vs `1.026`).
- **V16. One-sided limits at rationals (bounded-`p` sequences, `q ≤ 1153`, converged to 1e-5):**
  `x̃ → ⅓⁻` via `x̃=[0;3,N]` (`p=3`): `κ → 0.05453 + 0.00506i`, `G → 1.12837`.
  `x̃ → ⅓⁺` via `[0;2,1,N]` (`p=3`): `κ → 0.05504 − 0.02351i`, `G → 1.12537`.
  ⇒ at `x̃ = ⅓`: `Im κ` **jumps** by `−0.0286`, `Re κ` by `+0.0005`; `G` jumps by `−0.0030` (not a V-cusp alone: cusp + jump).
  `x̃ → ½⁻` via `[0;2,N]` (`p=2`): `κ → 0.04969 + 0.0240i`, `G → 1.11142` (supersedes the extrapolated `≈1.102` of V6/E-d).
  `x̃ → 0⁺` via `[0;N]` (`p=1`): `Re κ → 0.02383`, `Im κ → −0.0527 ± 0.001` (log-slow), `G → 1.0263–1.0265`.
- **V17. The limit depends on the far end of the continued fraction (order of limits).** Same `x̃` to `1e-5`, different `G`:
  `x̃ = [0;3,N,5]` (`p/q = [0;5,N,3]`, `t→⅕`): `G = 1.12551, 1.12462, 1.12371, 1.12326` at `N = 48, 64, 96, 128`, fit `G_∞ + 0.17/N` ⇒ `G_∞ ≈ 1.1219`, vs `1.12837` for `[0;3,N]`.
  `x̃ = [0;N,4]` (`p/q=[0;4,N]`, `t→¼`): `G → ≈1.014` (`1.01509` at `q=1537`, slope `0.5/N`) vs `1.0265` for `p=1`.
  At fixed `N`, `a₁ → ∞` (`p` fixed `= 3N+1`, `q → ∞`) converges smoothly: `[0;3,16,a₁]`: `G = 1.13089, 1.13164, 1.13179, …, 1.13143` for `a₁ = 2 … 100`; `[0;64,a₁]`: `1.02585 → 1.01689`. Grid `[0;3,N,a₁]`, `N ∈ {32,64,128}`, `a₁ ≤ 64`: see `data/limit_grid.json` (GRID_SUMMARY).
- **V18. The second variable is arithmetic in `t`, not smooth.** At tail `(N=128, 3)` (`x̃ = ⅓ − 9e-4`), prefixes giving `t = 0.5, 0.4, 0.4286, 0.4545, 0.333, 0.286, 0.25, 0.2, 0.111`: `G = 1.1260, 1.1238, 1.1232, 1.1241, 1.1246, 1.1236, 1.1238, 1.1233, 1.1238` — peaks at `t = ½` (`+0.0024`), `⅓` (`+0.001`) over a `≈1.1236` floor, and `t→0` with `p=3`: `1.1284`. Symmetry `t ↔ 1−t` exact (`(2,2)≡(1,1,2)`, `(3,)≡(1,2)`, `(4,)≡(1,3)`). Dense `q=1009` additive model `G ≈ F(x*) + H(t)`: DENSE_SUMMARY.
- **V19. Universality across components is exact for the `x̃`-part and fails only for bounded `p`.** `disk2/main2` ratio of `G` on `[0;a₁,16,3]`: `0.9924` at `p=3` (`a₁=1`), `1.0000 ± 0.0002` for `a₁ = 2…100`; on `[0;a₁,64]`: `0.9881` at `p=1`, `0.9997–1.0008` for `a₁ ≥ 2`. Refines V9 (`0.9985 ± 0.0032` at `q=59`, whose scatter was the bounded-`p` bulbs).
- **V20. Quadratic irrational:** `x̃ → 2−φ = [0;2,1,1,…]`: `κ = 0.05716, 0.06216, 0.06390` at `q = 34, 89, 610`; geometric convergence (ratio `≈ 1/3` per two Fibonacci steps), `κ_∞ ≈ 0.0648 − 0.0058i`, `G_∞ ≈ 1.158` (the global maximum region).

Precision caveats: doubles for `R_q`; the FFT route is accurate to `~1e-12` in `r₂` (verified against balls, V12) and is the instrument for `q > 300`. `G` from `|u_a|/2` carries the `O(q⁻²)` proxy difference of V14 (irrelevant for limits). Prime `q` used in sweeps; sequences use whatever `q` the CF gives.

---

## 4. FALSIFIED (this branch)

- F1 (was C12). `c₁(t)/q` correction with smooth `c₁`. **False** (V2–V3).
- F2 (was C13). "Nothing arithmetic transfers." **False** (V3–V6, V15).
- F3 (was C15). Same `Ĝ` for every degree. **False** (V10).
- **F4 (was C1″/C14′ as stated).** "`d = 2|φ'|q⁻² Ĝ(x*)(1+o(1))` with a single continuous function `Ĝ` on `(0,½]`." **False at the 0.5–1.2 % level**: the `q→∞` limit along `x̃ → x` depends on the sequence (V17), one-sided limits at rationals differ (V16: jump `0.003` at `⅓`). The leading `x*`-description survives as the dominant term (range 0.12 vs residual `≤ 0.013`).
- **F5 (was the contraction form of C16).** "Deeper ancestors enter with weight `(q_{n−1}/q)²`." **False**: the digit `a₁` of `p/q` (deepest ancestor of `x̃`) changes the limit by `0.005–0.012` even as the intermediate digit `N → ∞` (V17), which would be suppressed by `(1/N)²` under C16.
- ~~V6's claim of equal one-sided limits at ⅓~~ (superseded; the `0.001` agreement there was at one fixed prefix and inside its error).

---

## 5. CONJECTURED (open register)

**C1‴ (size law, corrected form).** For every hyperbolic component `H` of the quadratic family,
`d_H(p/q) = 2|c_H'(λ₀)| q⁻² · G_H(p/q) · (1 + O(q⁻²))`, with
`G_H(p/q) = Ĝ(x̃) + Ĥ_H(t; x̃) + o(1)`, where `Ĝ` is universal (component-independent) and depends on `x̃ = p̄/q` only, range `[1.014, 1.16]`; `Ĥ_H` is `≤ 0.013` in size, component-dependent, supported on bulbs with small-denominator `t = p/q` (bounded `p` in the limit), and `o(1)` is `O(1/N)` in the last partial quotient. Evidence: V15–V19.

**C14″ (structure of Ĝ and of κ̂).** `κ̂(x̃) = lim κ` (along sequences with `p → ∞`) exists for every irrational `x̃` and every one-sided approach to a rational; `Re κ̂` is continuous on irrationals with V-cusps at rationals; `Im κ̂` is odd, with **jumps** at rationals `p'/q'`; both singular amplitudes decay in `q'`. Proposed Fourier form `κ̂(x̃) = Σ_{m≥1} W_m e^{2πimx̃}/m^s + c.c.`-type with `s ∈ {1,2}` mixing Bernoulli (`B₂`, cusps) and sawtooth/Clausen (jumps) parts — the `m = kq'` terms produce the Ford weights. Open: exact cusp/jump laws in `q'` (V7 gave `≈ 0.18/q'²` for `Ĝ` on four points).

**C16′ (two-ended resonance structure, replaces C16).** `ι_{p/q}/q → Σ_{j} w(j)/(1 − λ₀^j)`-type sum over near-resonances `‖jp/q‖`. Two families of `j` are `O(1/q)`-resonant: `j ≡ mp̄ (mod q)` (`‖jp/q‖ = m/q`; gives `Ĝ(x̃)`, universal, "renormalised") and `j = k` small when `p` is bounded (`‖kp/q‖ = kp/q`; gives `Ĥ`, component-dependent because the first iterates see the component's own geometry). For `p = 1` both families coincide (explains why `p=1` is the extreme case: `G(1/q) → 1.0265` vs `≈1.017` for `[0;a₁,N]`, `a₁ → ∞`).

**C17′ (analytic mechanism, now half proven).** ~~`ρ = 1 − q²ε + κ q⁴ε² + …` with `κ` a small-denominator sum~~ → **P2**: `κ = (ι − ½)/q`. Remaining: the asymptotics `ι_{p/q} = ½ + q κ̂(x̃) + O(1)` as `q → ∞`, `p → ∞`. Route: index formula on the `W = w^q` quotient map `Φ(W) = W h(W)^q` (`ι = [W²]`-data of `Φ`), with `h` from the resonant normal form of `f_{λ₀}` whose coefficients are products of `1/(1 − λ₀^j)`; equivalently, `ι = −Σ_{other cycles of period | q} 1/(1−ρ)` [Mil] — a sum over the repelling cycles at the root, dominated by the least repelling ones. Cross-reference [BE02] for an a-priori bound `|Re ι| ≲ q`.

**C20 (`p=1` constant).** `κ₀ := lim_{q→∞} (ι_{1/q} − ½)/q = 0.0238 − 0.0527i` (Im to `±0.001`) is an Écalle–Voronin / horn-map invariant of `z + z²` (parabolic implosion at `λ₀ → 1`: the `1/q`-cycle's multiplier is governed by the Lavaurs phase), hence expressible via the Fatou coordinate of `z+z²`. Untested.

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

1. **C17′ → asymptotics of `ι`.** Compute `ι` through the cycle-index sum `−Σ 1/(1−ρ)` over the repelling cycles of `f_{λ₀}` (all `2^q` fixed points of `f^q` for `q ≤ 20`, exact in `ℚ(ζ_q)` via `scripts/exact_index.py` + dynatomic factorisation); identify the dominant cycles (expect: the cycles of the neighbouring bulbs `B_{p'/q'}` with `|pq' − p'q| = 1`, i.e. the Farey neighbours — which would be the Ford-circle mechanism made literal).
2. **Cusp/jump laws.** From `data/kappa_q1009.json`, fit `Re κ̂` cusps and `Im κ̂` jumps at `p'/q'` for `q' ≤ 7` against `q'^{-1}`, `q'^{-2}`; confirm with `q = 2003`.
3. **Fatou-coordinate computation of `κ₀`** (C20), cheap via the standard `z + z²` Fatou coordinate series.
4. **`z³+c`**: repeat V15–V19 (`MAIN3` is already wired: `G(MAIN3, p, q)`, `taylor(p, q, MAIN3)`).
5. C2′ numerics, then C8, C11 as before.

## 8. Instruments (`bulbford/`, tests in `tests/`, 17 passing)

- `cf.py` — `modinv`, `xstar`, `cf`, `from_cf`, `convergent_denominators`.
- `dynamics.py` — `Family` records (`MAIN2`, `DISK2`, `MAIN3`), `orbit` with second derivatives, `Cycle` tracking with analytic `dρ/dc` and a period-collapse guard, `bulb(fam, p, q)` → `G_ant`, `G_cen`; `rho_on_path`.
- `index.py` — `index(p, q)`: `ι_{p/q}` as an `acb` ball (python-flint), auto precision; `kappa(p, q)`.
- `taylor.py` — `taylor(p, q, fam, r, N)`: `r_k`, `solve(target, u0)`; `kappa_fft` (N=64) for `q > 300`.
- `scripts/` — `exact_index.py` (ℚ(ζ_q)), `taylor_sweep.py`, `index_sweep.py`, `kappa_sweep.py` (dense, `q=1009`), `sequences.py`, `prefix_test.py`, `a1_scan*.py`, `limit_grid.py`, `analyze.py`, `analyze_dense.py`. Data in `data/` (JSON). Legacy instruments in `legacy/`.
- Run: `PYTHONPATH=. pytest -q`; `PYTHONPATH=. python3 scripts/analyze.py 59 127 251`.

## 9. Error / caveat catalog

- E-a. Diameter proxy (ρ=−1 point). Unchanged; `|u_a|/2` is the cleaner proxy (V14).
- E-b. V7 rests on four cusp depths; the jump/cusp laws are open (next move 2).
- E-c. `q > 2100` untested; period-collapse guard bisects continuation steps (no failures logged).
- E-d. ~~`Ĝ(½)` extrapolated~~ → direct limit `1.11142` along `p=2`; other approaches to `½` differ at the V17 level.
- E-e. "Converged" for sequences means stable to `1e-5` across a doubling of `N`; the `[0;…,N,a]` families converge only like `c/N` and are quoted with the fit.
- E-f. The `ι`-route needs `dps ≈ 1.5q` and `O(q³)` ball operations; used to `q = 251` (385 s). Beyond that the FFT route is the instrument (validated to `1e-12` by V12).
