# RESEARCH BRANCH — Bulb geometry of M vs Ford circles: the arithmetic correction Ĝ

Branch id: `bulb-ford-correction`  ·  opened 2026-09-24  ·  status: **active, numerics-validated, no theorems yet**
Discipline: PROVEN / VALIDATED (numerical, stated precision) / CONJECTURED / FALSIFIED are absolute. Analytic imports carry theorem tags.

---

## 0. One-paragraph summary

The satellite bulbs of the main cardioid (and of every hyperbolic component of the quadratic family) obey
`d(p/q) = 2|φ_H'(e^{2πip/q})| · q⁻² · Ĝ(‖p⁻¹ mod q‖ / q) · (1+o(1))`
with a single universal function `Ĝ : (0, 1/2] → [1.026, 1.15]`, continuous, with a downward V-cusp at every rational `p'/q'` of depth `≍ q'⁻²`. The Ford-circle analogy is exact at the combinatorial level (Farey adjacency) and at leading order (`q⁻²`), and at first correction the Ford radius law reappears as the *singular structure* of `Ĝ`. The modular group enters at first correction arithmetically (`p ↦ p⁻¹ mod q`), not geometrically. Universality holds across hyperbolic components of `z²+c` but **fails across degree** (`z³+c` has its own `Ĝ₃`).

---

## 1. Definitions & notation

- `B_{p/q}` : hyperbolic component of period `q` attached to the main cardioid at internal angle `p/q`; root `c_root(p/q) = λ₀/2 − λ₀²/4`, `λ₀ = e^{2πip/q}`.
- `d(p/q)` : diameter proxy `|c_ant − c_root|`, where `c_ant` is the boundary point of `B_{p/q}` with satellite multiplier `ρ = −1` (root of the period-2q satellite), reached by continuation from the center along the ray root→center. Secondary proxy `2|c_cen − c_root|`.
- `φ_H` : inverse of the multiplier map `λ_H : H → 𝔻` of a hyperbolic component `H`. Cardioid: `φ' = (1−λ)/2`, `|φ'(λ₀)| = sin(πp/q)`. Period-2 disk: `φ = λ/4 − 1`, `|φ'| = 1/4`.
- `x*(p/q) := q_{n−1}/q = [0; a_n, …, a_1]` (reversed continued fraction) `= ‖p⁻¹ mod q‖/q = min(p̄, q−p̄)/q`, `p·p̄ ≡ 1 (mod q)`. Ranges over `(0, 1/2]` for canonical CF with `a_n ≥ 2`.
- `G(p/q) := d(p/q) / (2|φ'(λ₀)| q⁻²)`; `Ĝ` its limit as a function of `x*`.
- Ford circle `C(p/q)`: tangent to ℝ at `p/q`, radius `1/(2q²)`; `C(p/q) ⟂ C(r/s)` tangent iff `|ps−qr| = 1`.

Theorem tags used: **[GM84]** Guckenheimer–McGehee `diam B_{p/q} = O(q⁻²)`; **[DH]** Douady–Hubbard multiplier map / external angles of roots; **[BS94]** Bullett–Sentenac, external angles of `B_{p/q}` roots are Christoffel-word rationals; **[Yoc95]** Yoccoz `log r(α) = −B(α) + O(1)`; **[MMY]** Marmi–Moussa–Yoccoz ½-Hölder conjecture; **[Weil]** Kloosterman bound / equidistribution of `(p/q, p̄/q)`; **[FL]** Franel–Landau.

---

## 2. PROVEN (this branch)

Nothing beyond the trivial:
- P0. First-order local model: with `f^q(z) = λ^q z + a z^{q+1} + …`, the period-q cycle multiplier is `ρ = 1 − q(λ^q − 1) + O((λ^q−1)²)`; for `λ = λ₀e^ε`, `ρ = 1 − q²ε + O(q⁴ε²)`. Hence the `ε`-disk of radius `q⁻²` and the leading law `d = 2|φ'(λ₀)| q⁻²` **to first order**. (Elementary; the `O(q⁴ε²)` term is *not* small at the boundary — see V2.)
- P1. `q = 2` exact: period-2 disk diameter `1/2 = 2·sin(π/2)/4`.

---

## 3. VALIDATED (numerical; doubles; Newton to 1e-14 relative; `q ≤ 2014`)

- V1. Leading constant `2` and factor `|φ'(λ₀)|` confirmed: `G ∈ [0.98, 1.03]` for `q ≤ 5`; `G → Ĝ ∈ [1.026, 1.15]` for large `q`.
- V2. **`G` does not tend to 1.** `q·(G−1)` grows linearly in `q`; the correction is `O(1)` multiplicative, not `O(1/q)`. The first-order disk model is *not* asymptotically exact in shape.
- V3. **`G` is not a function of `t = p/q`.** At `q = 59`, `G` jumps between 1.026 and 1.145 as `p` steps by 6. Correlations at `q = 127`: `|Dedekind s(p,q)|` −0.51, `Σaᵢ` −0.54, `max aᵢ` −0.54 — none deterministic.
- V4. **`G` is a near-deterministic function of the last partial quotient `a_n`** (`q ∈ [40,100]`, ~1000 bulbs, cf length ≥ 3):
  `a_n=2: 1.141(sd .009)  3: 1.135  4: 1.122  5: 1.109  6: 1.097  8: 1.080  10: 1.068  12: 1.060  15: 1.051  19: 1.044  →∞: 1.026`.
  Residual dependence on `a_{n−1}` ~0.005 (an order of magnitude smaller).
- V5. **`G` is continuous in `x*`.** Sorted by `x*` at `q = 59`: smooth unimodal curve, `1.026` at `x*→0`, max `≈1.147` near `x* ≈ 0.42`, `≈1.11` at `x* → 1/2`; max adjacent jump 0.015 vs range 0.12. Sample (x*, G_ant, G_cen): (0.0169, 1.0262, 1.0195) (0.1017, 1.0725, —) (0.2034, 1.1112, —) (0.3051, 1.1384, —) (0.4237, 1.1471, —) (0.4915, 1.1113, —).
- V6. **One-sided limits at rationals agree; V-cusp.** `x* → 1/3` via tails `(2,1,N)` and `(3,N)`, prefix `[1,2,1,2]`, `N ∈ {3,…,60}`, `q` up to 2014: both → `≈1.120` (fit `G∞ + c/N`: `G∞ = 1.1198`, `c ≈ 0.18` from above; `≈1.121` from below). Envelope near 1/3 ≈ 1.138 ⇒ cusp depth ≈ 0.018. At `x*→1/2`: → ≈1.102 (envelope ≈1.147, depth ≈0.045). At `x*→0`: → 1.026 with linear approach (slope ≈0.6).
- V7. **Cusp depths vs denominator** (from V5/V6, q=59 and limits): `1/2: 0.045, 1/3: 0.018, 1/4: 0.010, 1/5: 0.005–0.007` — consistent with `depth ≈ 0.18/q'²` to within noise. *Four data points; weakest validated item.*
- V8. **Asymmetry** of the bulb (`|c_ant−c_root| / 2|c_cen−c_root|`) tracks `x*`: 1.007 (`x*→0`) to 1.078 (`x*≈0.42`). The limit shape is a function of `x*`.
- V9. **Universality across components** (`q=59`): period-2-disk satellites vs cardioid, ratio `G_disk2/G_main2 = 0.9985 ± 0.0032`, corr 0.996. Same `Ĝ`.
- V10. **Non-universality across degree**: `z³+c` main component, corr with `Ĝ₂` 0.95 but ratio `1.081 ± 0.017` and `Ĝ₃(0⁺) ≈ 1.171`, peak `≈1.22`. Different function, not a scalar multiple.
- V11. Sum of diameters: `Σ_{q≤Q} d ≈ 0.83–0.89 · log Q` vs first-order prediction `24/π³ ≈ 0.774`; ratio consistent with `⟨Ĝ⟩ ≈ 1.08–1.09`.

Precision caveats: doubles; `d` down to ~3e-7 at `q≈2000` with absolute Newton tolerance ~1e-14·d — adequate for 3–4 significant digits in `G`. `c_ant` is a diameter *proxy* (ρ=−1 point, not sup-diameter). Prime `q` used to avoid degenerate `p`.

---

## 4. FALSIFIED (this branch)

- F1 (was C12). `d(p/q) = 2 sin(πp/q)q⁻²(1 + c₁(t)/q + O(q⁻²))` with `c₁` smooth in `t`. **False** by V2–V3.
- F2 (was C13, the "nothing arithmetic transfers" meta-conjecture). **False** by V3–V6: second-order bulb geometry is a continuous function of the modular inverse.
- F3 (was C15 as stated, "same `Ĝ` for every degree"). **False** by V10; survives only within a family (V9).

---

## 5. CONJECTURED (open register)

**C1″ (size law, final form).** For every hyperbolic component `H` of `M`:
`d_H(p/q) = 2|φ_H'(e^{2πip/q})| q⁻² Ĝ(x*) (1 + o(1))`, `Ĝ` universal for the quadratic family. Evidence: V1, V4–V6, V9.

**C14′ (structure of Ĝ).** `Ĝ` is continuous on `[0,1/2]`, smooth off ℚ, with a V-cusp at every rational `p'/q'` of depth `≍ q'⁻²`, one-sided Lipschitz at each cusp. Proposed form `Ĝ(x) = Ĝ₀(x) − Σ_{p'/q'} q'⁻² V(q'²(x − p'/q')) + …` (Takagi/Brjuno-type). Open: Hölder exponent of `Ĝ`; whether `Ĝ` is differentiable anywhere on ℚ-closure. Evidence: V5–V7.

**C16 (renormalization functional equation).** `log Ĝ(x) = U₁(x) + x² U(Tx) + …` with `T` the Gauss map, i.e. `Ĝ` is a fixed point of a transfer-type operator with contraction `(q_{n−1}/q)²` per ancestor. Explains geometric tail-decay (V4 residuals) and Ford-weighted cusps (V7).

**C17 (analytic mechanism).** Expand the satellite multiplier `ρ = 1 − q²ε + κ(p/q) q⁴ε² + …`. Then `Ĝ` is an explicit function of `κ`, and `κ` is a small-denominator sum of type `Σ_{j<q} w(j)/(1−e^{2πijp/q})²`, dominated by near-resonant `j ∈ q_{n−1}ℤ` (`‖jp/q‖ = m/q`). Continuity in `q_{n−1}/q` and cusps where the next convergent takes over follow; Dedekind sums (uniform weights) do not see this. **This is the route to proving C1″/C14′.**

**C2′ (bulb zeta).** `Z_M(s) := Σ_{p/q} d(p/q)^s = 2^s ζ(2s−1)/ζ(2s) · ∫₀¹ sin^s(πt)dt · ⟨Ĝ^s⟩ + H(s)`, `Res_{s=1} = (12/π³)⟨Ĝ⟩`, `⟨Ĝ⟩ ≈ 1.08–1.09`. Continuation of `H` past `Re s = 1` is governed by Kloosterman sums via [Weil] equidistribution of `(p/q, p̄/q)`. Evidence: V11.

**C9′ (limit measure).** `μ_Q := Σ_{q≤Q} Σ_p d(p/q) δ_{p/q}`, `μ_Q / log Q → (12/π³)⟨Ĝ⟩ sin(πt) dt`; rate is a Kloosterman statement.

**C8 (Christoffel discrepancy).** With `θ_−(p/q)` the external angles of roots [BS94], the discrepancy of `{θ_−(p/q)}_{q≤Q}` has exponent strictly in `(0, 1/2)` and is not RH-sensitive (the [FL] mechanism needs the identity map on Farey fractions). Untested.

**C11 (quadratic irrationals).** For `α` quadratic irrational with primitive hyperbolic `γ_α`, eigenvalue `ε_α`: rescaled bulb shapes along convergents converge with rate `O(ε_α^{−n})`; the `O(1)` term in [Yoc95] restricted to quadratic irrationals is a function of the CF period. Untested; connects to [MMY].

**C15′ (family-universality, corrected).** `Ĝ` is universal across all hyperbolic components of one family (V9) and is a functional of the parabolic-germ data of the family beyond the quadratic term (V10). Open: is the cusp law `q'⁻²` family-universal even though the envelope is not?

---

## 6. Relationship to Ford circles (current statement)

| level | Ford circles | bulbs |
|---|---|---|
| combinatorics | Farey adjacency, mediants | identical |
| leading order | radius `1/(2q²)`, exact | `2|φ'|q⁻²`, asymptotic |
| first correction | none (PSL(2,ℤ)-exact) | `Ĝ(‖p⁻¹/q‖)`; singular set ℚ with Ford weights `q'⁻²` |
| symmetry | PSL(2,ℤ) on horoballs | PSL(2,ℤ) via `p ↦ p⁻¹` on the boundary |

---

## 7. Next moves (ordered)

1. **C17 → C1″.** Compute `κ(p/q)` exactly from the resonant normal form of `λz+z²` at `λ = e^{2πip/q}`; verify against V5 table; obtain `Ĝ` in closed form. This is the tractable theorem.
2. **Hölder exponent of `Ĝ`.** Dense `x*` sampling at `q ≈ 500–1000` (prime), local oscillation vs scale. Cheap.
3. **C15′.** Recompute the cusp law for `z³+c`: if `q'⁻²` persists with a different envelope, the Ford weights are family-universal.
4. **C2′ numerics.** `Σ d` to `Q ≈ 500` with the fitted `⟨Ĝ⟩`; test the `ζ(2s−1)/ζ(2s)` structure via Möbius inversion on `q`.
5. C8 (pure combinatorics on Christoffel words), then C11.

## 8. Instruments

- `bulbs.py` — cardioid bulbs: center (Newton on `f^q(0)=0`), antipode (continuation + Newton on `ρ=−1`), tables, `Σd` vs `24/π³ log Q`.
- `general.py` — family-generic (`main2`, `disk2`, `main3`): `G(family, p, q)` → `(G_ant, G_cen)`.
- Run `python3 bulbs.py`; `python3 general.py` (q=59 cross-family table). `from general import G` for ad-hoc CF-constructed `p/q`.

## 9. Error / caveat catalog

- E-a. Diameter proxy is the ρ=−1 point; sup-diameter may differ by an `x*`-dependent factor (bounded by V8 asymmetry ≈ 8%). Does not affect any qualitative claim; affects `Ĝ` values at the 1e-2 level.
- E-b. V7 rests on four cusp depths; treat `q'⁻²` as a hypothesis, not a fit.
- E-c. Doubles; `q > 2000` untested. Newton continuation could in principle track a wrong cycle branch — guarded by checking `|ρ+1| < 1e-8` (no warnings fired).
- E-d. `Ĝ(1/2)` is a one-sided limit (`a_n ≥ 2`); the endpoint value ≈1.102 is extrapolated.
