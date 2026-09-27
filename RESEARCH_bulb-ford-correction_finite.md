# RESEARCH BRANCH, finite version: bulb geometry of M vs Ford circles, without limits

Branch id: `bulb-ford-correction`  ·  finite revision 2026-09-27 of `RESEARCH_bulb-ford-correction.md` (2026-09-25)  ·  status: **active**
Discipline: PROVEN / VALIDATED / CONJECTURED / FALSIFIED are absolute, as in the original register. This version adds a second rule. **No claim below is a limit.** Each claim of the original takes one of the statement forms of §0, or moves to §5 as a named classical referent that this document never asserts. The original register is unchanged and remains the claim-state file of record. §8 maps every original item to its form here.

---

## 0. Statement forms

| Form | What it asserts | Replaces |
|---|---|---|
| **E** exact identity | an equation in a named ring: `ℚ(ζ_q)`, `ℚ`, or the convergent power series `ℂ{ε}` | Taylor remainders, "total index converges" |
| **C** certificate | a finite check in rational interval or ball arithmetic, replayable from stored data | "the root is near" |
| **T** finite table | values at stated parameters, with the stated precision of the instrument | "`→ L` as `q → ∞`", "`O(q⁻ᵏ)`" |
| **S** fitted statistic | the least-squares parameters of a stated model on stated nodes | extrapolated limits (`G_∞`, `κ_∞`) |
| **X** finite falsification | a finite set of data points that no member of a stated function class fits | "`G` is not a function of `x̃` alone" |
| **K** classical referent | a named statement about infinitely many bulbs, recorded with its finite shadow and never asserted | the conjectures C1‴ … C20 |

Two replacements carry the load.

1. **Coefficient extraction replaces the remainder.** `f(ε) = a + bε + O(ε²)` becomes `[ε⁰]f = a`, `[ε¹]f = b` for a convergent power series `f ∈ ℂ{ε}`. A Laurent coefficient is read off an identity of series, not approached.
2. **A point of `ℙ¹` replaces the value at infinity.** A value at `x = ∞` is the value at `[1:0]` in the chart `s = 1/x` (`larsbx/finite-math-kernels`, `projective_limits`). An S-statistic is a named finite computation, and its intercept is a fit parameter, never a limit.

**Enforcement.** `tools/audit_limits.py` (run in CI and by `tests/test_audit_limits.py`) rejects the idioms `→`, `lim`, `O(`, `o(`, "limit", "converges", "asymptotic", "tends to" and "approaches" everywhere in this file except §0, §5 and §8, a clause that denies the idiom, or a paragraph marked `<!-- limit-exempt: reason -->`. Continued-fraction convergents and convergent power series are not limit idioms.

---

## 1. Definitions & notation

As in the original §1, with these changes.

- `G(p/q) := d(p/q) / (2|c_H'(λ₀)| q⁻²)` is a number attached to one bulb. No law `d = … · (1 + O(q⁻²))` is asserted; `G` is the finite ratio itself.
- `r_k(p/q) := [u^k] R_q(u)`, `u = q²ε`. `κ(p/q) := (ι_{p/q} − ½)/q` exactly.
- `u_a`: the root of `R_q(u) = −1` found by Newton from `u = 2`. `|u_a|/2` is a second proxy for `G`, compared to `G_ant` in V14 as a table only.
- `x̃ := p̄/q ∈ (−½, ½]`, `x* := |x̃|`, `t := p/q`, all exact rationals (`bulbford/cf.py`; R1 vectors of `finite-math-kernels`).
- Tags as in the original, plus **[Hol]**: a contour integral of a holomorphic family is holomorphic in the parameter; a holomorphic function is determined by its Laurent expansion. **[GM84]** is recorded but used by no claim here.

---

## 2. PROVEN

- **P0 + P2 (one Laurent identity), form E.** For any unicritical family, any hyperbolic component `H` and `λ = λ₀ e^ε`, with `μ = λ^q = e^{qε}` and `ρ(ε)` the satellite multiplier, `ρ ∈ ℂ{ε}` satisfies
  `[ε⁰]ρ = 1`,  `[ε¹]ρ = −q²`,  `[ε²]ρ = q³(ι − ½)`,
  where `ι = ι(F, z₀)` is the holomorphic index of the parabolic point of `F = f_{λ₀}^{q·per(H)}`.
  *Proof.* Fix a circle `γ` around `z₀` that encloses no other fixed point of `F₀`. Put `I(ε) := (1/2πi)∮_γ dz/(z − F_ε(z))`. By [Hol], `I ∈ ℂ{ε}`, and by definition of the index `[ε⁰]I = I(0) = ι`. `ρ ∈ ℂ{ε}` as in the original (a symmetric function of the cycle points). The number of fixed points inside `γ` is the argument-principle integral, an integer-valued element of `ℂ{ε}`, hence the constant `q + 1`. By [DH] they are the `H`-point (multiplier `μ`) and the `q` satellite points (multiplier `ρ`). On a punctured neighbourhood of `ε = 0` all are simple, since `μ ≠ 1` and `ρ ≢ 1` has isolated zeros. The residue theorem gives the identity of Laurent series
  `1/(1−μ) + q/(1−ρ) = I(ε)`.
  Here `1/(1−μ) = −1/(qε) + ½ − qε/12 + …`. Write `ρ = 1 + ρ_k ε^k + …` with `ρ_k ≠ 0`. The right side has no pole, so the left side has none. For `k ≥ 2`, the pole of order `k` in `q/(1−ρ)` has nothing to cancel it, so `k = 1`. Cancelling the simple pole, `−1/q − q/ρ₁ = 0`, gives `ρ₁ = −q²`. The constant term gives `½ + q ρ₂/ρ₁² = ι`, that is `ρ₂ = q³(ι − ½)`. ∎
  Nothing is approached. The proof compares coefficients of two equal elements of `ℂ((ε))`, and "the total index of the merging points" is the constant term `[ε⁰]I`.
  *Exact checks.* `q = 1`: `ρ = 2 − e^ε`, so `[ε²]ρ = −½ = 1·(0 − ½)`, with `ι(z+z²) = 0`. `q = 2`: `ρ = 1 − 4ε − 3ε² + …`, `−3 = 8(⅛ − ½)`, with `ι = ⅛`.
  *Remarks.* The remarks of the original P2 stand. Remark (ii), `ι_{p'/q} = σ(ι_{p/q})` in `ℚ(ζ_q)`, is a form-E statement.
- **P1.** `q = 2`: period-2 diameter `1/2`, and `ρ(λ) = −λ² + 2λ + 4` for the 2-cycle of `λz + z²` (form E, tested).
- **P3 (certified centres), form C.** As in the original: 79 boxes, `2 ≤ q ≤ 16`, `data/center_certificates.json`, replayed by `tests/test_certify.py`.
- **P4 (wake combinatorics), form E.** As in the original: exact `ℚ/ℤ` arcs, Farey order, mediant in the gap, for `q ≤ 16`.
- **P5 (exact `ι`, `q ≤ 8`), form E.** As in the original: cyclotomic vectors replayed coordinate by coordinate, and `ι_{p/q} = σ_p(ι_{1/q})` holds exactly.
- **P6 (projective index sum), form E.** For a polynomial `F` of degree `d ≥ 2`, the lift `[Z:W] ↦ [W^d F(Z/W) : W^d]` fixes `∞ = [1:0]` with multiplier `0`, so `ι_∞ = 1`. The index formula on `ℙ¹` [Mil, §12], `Σ_{z ∈ ℙ¹} ι_z = 1`, gives `Σ_{z ∈ ℂ} ι_z = 0`. Apply it to `F = f_{λ₀}^{q}`. Every fixed point `z ≠ z₀` lies on a cycle of period dividing `q`, and these cycles are repelling [DH: a parabolic quadratic map has one non-repelling cycle]. So each such `z` is simple with `ι_z = 1/(1 − ρ_z)`, and
  `ι_{p/q} = −Σ_{fixed points z ≠ z₀} 1/(1 − ρ_z)`,
  as an exact identity at `ε = 0`, over the finitely many other fixed points. This is the finite form of the route in C17′: it evaluates `ι` as a finite sum rather than as the index of a merger.
- **P7 (flank-orbit lemma), form E.** Let `0 < p < q` with `gcd(p, q) = 1` and `q ≥ 3`. Put `M = 2^q − 1` and `p̄ = p⁻¹ mod q`. Let `α_0 < … < α_{q−1}` be the numerators over `M` of the `p/q` rotation cycle, and `(α_w, α_{w+1})` its characteristic arc. Then `w = p − 1` and, with indices mod `q`,
  `α_{p+1} − 1 ≡ 2^{p̄} (α_{p−2} + 1)  (mod M)`,
  and `α_{p−2} + 1` has exact period `q` under doubling. So the two flank angles of V22 lie on one doubling orbit of period `q`, and the second is reached from the first after `p̄` doublings.
  *Proof.* (i) *Words.* A numerator mod `M` is a cyclic `q`-bit word `W` (bit `0` most significant), and doubling is the left rotation `(σW)_k = W_{k+1}`. For `r ∈ ℤ/q` put `c(r)_k = [(r + kp) mod q ≥ q − p]`. Then `c(0)` is the word `rotation_cycle` starts from and `σc(r) = c(r + p)`, so the cycle is `{c(r)}`, and `σ^{p̄} c(r) = c(r + 1)`. (ii) *Neighbours.* `s_k = (r + kp) mod q` is a bijection of `ℤ/q`. Passing from `r` to `r + 1` changes bit `k` exactly when `s_k = q − p − 1` (from `0` to `1`) or `s_k = q − 1` (from `1` to `0`). The first happens at `k₁ ≡ (q − p − 1 − r) p̄`, and then `s_{k₁+1} = q − 1`. So `c(r+1)` is `c(r)` with its pair at `(k₁, k₁+1)` changed from `01` to `10`. (iii) *Order.* For `0 ≤ r ≤ q − 2`, `k₁ = q − 1` would force `r ≡ −1`, so the pair does not wrap, and `N(c(r+1)) − N(c(r)) = 2^{q−2−k₁} > 0`. Hence `α_r = N(c(r))`. The difference is `1` iff `k₁ = q − 2` iff `r = p − 1`, so `w = p − 1`. (iv) *`2 ≤ p ≤ q − 2`, so `p̄ ∉ {1, q−1}`.* By (ii), passing from `c(p−1)` to `c(p)` flips the pair `(q−2, q−1)`, from `c(p−2)` to `c(p−1)` flips `(p̄−2, p̄−1)`, and from `c(p)` to `c(p+1)` flips `(q−p̄−2, q−p̄−1)`. Since `p̄ ∉ {1, q−1}`, the last two pairs neither wrap nor meet `(q−2, q−1)`. So `c(p−2)` ends in `01`, and `X := α_{p−2} + 1` is `c(p−2)` with its last pair changed to `10`, that is, `c(p)` with its pair at `(p̄−2, p̄−1)` changed from `10` to `01`. Likewise `c(p+1)` ends in `10`, and `Y := α_{p+1} − 1` is `c(p+1)` with its last pair changed to `01`. `σ^{p̄}` sends `c(p)` to `c(p+1)` and bit position `P` to `P − p̄`, which sends `(p̄−2, p̄−1)` to `(q−2, q−1)`. So `σ^{p̄} X = Y`. (v) *`p = 1`.* `c(r)` has a single `1`, at position `q − 1 − r`, so `α_r = 2^r`, `w = 0`, `X = 2^{q−1} + 1` and `Y = 3 = 2X mod M`, with `p̄ = 1`. *From `p` to `q − p`.* Conjugation `x ↦ M − x` commutes with doubling and reverses the order. It sends the `(q−p)/q` cycle to the `p/q` cycle with `α'_i = M − α_{q−1−i}`, hence `X' = M − Y` and `Y' = M − X`. So the identity for `p` gives it for `q − p`, with `q − p̄ = (q−p)⁻¹ mod q`. This covers `p = q − 1`. (vi) *Period.* A word of period `d | q`, `d < q`, has weight divisible by `q/d`. For `2 ≤ p ≤ q − 2` the weight of `X` is `p` (flips preserve weight), which is prime to `q`, so `d = q`. For `p = 1`, the two `1`s of `X` are cyclically adjacent, which is impossible for a period `d < q` when `q ≥ 3`. `p = q − 1` follows by conjugation. ∎
  *Checked in integers* (`tests/test_flank_lemma.py`, `bulbford.wake.mechanical`): each step (i)–(v) and the conclusion for all 1100 fractions with `q ≤ 60`, and agreement with `rotation_cycle`, `wake` and `flank_angles` for `q ≤ 16`.
- **P8 (second-ring lemma), form E.** In the setting of P7 with `q ≥ 5`, put `X₂ = α_{w−2} + 1` and `Y₂ = α_{w+3} − 1` (numerators over `M`, indices mod `q`, `w = p − 1`). Then `X₂` and `Y₂` lie on **distinct** doubling orbits, each of period exactly `q`. Unlike ring 1 (P7), ring 2 splits into two cycles.
  *Proof.* (i) *Height set.* For a cyclic `q`-bit word `W` of weight `m`, let `H_k` be the number of `1`s among its first `k` bits and `g_k = qH_k − km` (`0 ≤ k < q`). Extending `H` by `H_{k+q} = H_k + m` makes `g` periodic, and the rotation `σ^t W` has `g'_k = g_{k+t} − g_t`. So the set `G(W) = {g_k}`, taken up to translation, is constant on each doubling orbit. (ii) *Mechanical words.* For `c(r)`, `0 ≤ r < q` (P7), `H_k = ⌊(r + kp)/q⌋`, so `g_k = r − s_k` with `s_k = (r + kp) mod q` a bijection. Hence `G(c(r))` is the interval `[r − q + 1, r]`. (iii) *`3 ≤ p ≤ q − 3`.* Here `X₂ = N(c(p−3)) + 1` and `Y₂ = N(c(p+2)) − 1`, with `0 ≤ p − 3` and `p + 2 ≤ q − 1`. The last two bits of `c(p−3)` are `[s_{q−2} ≥ q−p] = [q−p−3 ≥ q−p] = 0` and `[s_{q−1} ≥ q−p] = [q−3 ≥ q−p] = 1`. So adding `1` changes only that last pair, from `01` to `10`. This raises `H_{q−1}` by `1`, so it raises `g_{q−1} = (p−3) − (q−3) = p − q` to `p`, and leaves every other `g_k` unchanged. Likewise `c(p+2)` ends in `10` (`s_{q−2} = q + 2 − p ≥ q − p`, `s_{q−1} = 2 < q − p`), and subtracting `1` lowers `g_{q−1} = p` to `p − q`. So `G(X₂) = {p−q−2, p−q−1} ∪ [p−q+1, p−3] ∪ {p}` and `G(Y₂) = {p−q} ∪ [p−q+3, p−1] ∪ {p+1, p+2}`. Both sets have `max − min = q + 2`, so a translation carrying one onto the other must match maxima, that is, shift `G(Y₂)` by `−2`. But `p − q − 1 ∈ G(X₂)` is not in `G(Y₂) − 2 = {p−q−2} ∪ [p−q+1, p−3] ∪ {p−1, p}`. So the orbits differ. Both words keep weight `p`, prime to `q`, so both have period `q` (P7 (vi)). (iv) *`p = 1`.* `X₂ = 2^{q−2} + 1` has weight `2` and `Y₂ = 7` has weight `3`, and weight is a doubling invariant. Their periods are `q`: two `1`s at cyclic distance `2` would force `q = 4`, and three consecutive `1`s would force `q = 3`. *`p = 2`* (`q` odd). `c(q−1)` has its `1`s at `0` and `(q−1)/2` and ends in `00`, so `X₂` has `1`s at `q−1, 0, (q−1)/2`, with cyclic gaps `(1, (q−1)/2, (q−1)/2)`. `c(4)` has its `1`s at `(q−5)/2` and `q−3` and ends in `100`, so `Y₂` has `1`s at `(q−5)/2, q−2, q−1`, with gaps `(1, (q−3)/2, (q+1)/2)`. The multiset of cyclic gaps is a rotation invariant, and these differ. Neither is constant, so both periods are `q`. *From `p` to `q − p`.* Conjugation `x ↦ M − x` commutes with doubling and exchanges the ring-2 angles of `p/q` and `(q−p)/q` (`X₂' = M − Y₂`, `Y₂' = M − X₂`). This covers `p = q − 2, q − 1`. ∎
  *Checked in integers* (`tests/test_ring_lemma.py`, `bulbford.wake.heights`): each step and the conclusion for all 1096 fractions with `5 ≤ q ≤ 60`. Using `m + 1` in place of `m` in `heights` fails 3120 of them.
- **P9 (arc-width lemma), form E.** In the setting of P7, the arc from `α_{w+j}` to `α_{w+j+1}` (indices mod `q`; for `w + j ≡ q − 1` this is the arc through angle `0`) has width `2^{e_j}/M` with `e_j = (j·p̄) mod q`. So Goldberg's `q` arc widths `2^0/M, …, 2^{q−1}/M` are placed around the characteristic arc by `x̃ = p̄/q`. In particular, the two arcs flanking the characteristic arc have widths `2^{p̄}/M` (`j = 1`) and `2^{q−p̄}/M` (`j = −1`). One of them is among the two widest (`e ≥ q − 2`, width at least `2^{q−2}/M`) if and only if `x* ≤ 2/q`.
  *Proof.* For `0 ≤ r ≤ q − 2`, P7 (iii) gives `α_{r+1} − α_r = 2^{q−2−k₁}` with `k₁ ≡ (q − p − 1 − r)·p̄` and `0 ≤ k₁ ≤ q − 2`. With `r = w + j = p − 1 + j`, `k₁ ≡ −2 − j·p̄`, so `q − 2 − k₁ ≡ j·p̄ (mod q)` and `0 ≤ q − 2 − k₁ ≤ q − 2`. Hence the exponent of each of these `q − 1` arcs is `e_j`, and `e_j ≠ q − 1` for them. Since `j ↦ j·p̄` is injective, their exponents are exactly `0, …, q − 2`. Their widths sum to `2^{q−1} − 1`, so the remaining arc, through `0`, has width `M − (2^{q−1} − 1) = 2^{q−1}`, the exponent `q − 1 = e_j` for its `j ≡ −p`. For the flanking arcs, `e_1 = p̄` and `e_{−1} = q − p̄` with `1 ≤ p̄ ≤ q − 1`. The larger is at least `q − 2` iff `p̄ ≤ 2` or `p̄ ≥ q − 2`, iff `x* = min(p̄, q − p̄)/q ≤ 2/q`. ∎
  *Checked in integers* (`tests/test_arc_widths.py`): every arc width and the flanking-arc equivalence for all fractions with `q ≤ 60`.
- **P10 (certified antipodes), form C.** For every `p/q`, `2 ≤ q ≤ 16` (79 cases), a box `Z × C` of half-width `2⁻⁶⁴` carries a two-variable Krawczyk inclusion for `(f_c^q(z) − z, (f_c^q)'(z) + 1)` and the type-`(0, q)` exclusions on the orbit of `Z`: exactly one `(z, c)` in the box, `z` of exact period `q` with multiplier `−1`. `ζ_q` is a Krawczyk box selected by an exact order check (no angle), `λ₀ = ζ_q^p` and `c_root = λ₀/2 − λ₀²/4` are interval expressions, and `G_ant² ∈ q⁴ Qd(C − c_root)/Qd(1 − λ₀)` is bracketed with width below `10⁻¹⁶` (`data/antipode_certificates.json`, replayed by `tests/test_antipode.py`). Reading `C` as the antipode of `B_{p/q}` uses [DH] and `SatelliteLabel` (§5 referents; the label is also checked as a table: at all 79 certified centres the critical orbit is cyclically ordered about `α` with step `p`, form T).

- **P11 (the β term in closed form), form E.** At the root of `B_{p/q}`, the fixed points of `z² + c` are `λ₀/2` (the parabolic point `z₀`) and `β = 1 − λ₀/2`: they sum to `1` and multiply to `c = λ₀/2 − λ₀²/4`. Their multipliers are `2z`, that is, `λ₀` and `2 − λ₀`. So `β`, the landing point of the ray of angle `0` [DH], contributes the exact term `T_β = −1/(1 − (2 − λ₀)^q)` to the P6 sum, with `|2 − λ₀|² = 5 − 4 cos(2πp/q) = 1 + 8 sin²(πp/q)`. ∎ *Checked* against the orbit computation of angle `0` (`tests/test_cycles.py`, 7 fractions up to `q = 64`, to `1e-9`).
---

## 3. VALIDATED (finite tables; doubles unless stated; Newton to `1e-14` relative)

- **V12 (P2 against instruments), form T.** The deviation `|[u²]R_q − κ|`, with the coefficient taken by Cauchy/FFT on `|u| = 2.5` (256 nodes) and `κ` in ball arithmetic of radius `< 1e-12`, has these maxima over all `p`:

  | `q` | 59 | 127 | 251 |
  |---|---|---|---|
  | max deviation | `9e-15` | `1e-13` | `1.3e-12` |

  `|r₀ − 1|, |r₁ + 1| < 1e-11` at the same `q`.
- **V13 (exact values), form E.** `ι_{1/3} = (92 − 16ζ₃)/441`, `ι_{1/4} = (1447 − 365i)/4624`, `ι_{1/5} = (33108 − 9153ζ − 7113ζ² − 312ζ³)/93775`. The first two denominators are squares of algebraic norms: `441 = 21²`, `4624 = 68²`. For `q = 3 … 12`, successive differences of `Im ι(1/q)` lie in `[−0.065, −0.063]` (`data/exact_index.txt`).
- **V14 (two proxies), form T.** `Δ := |(|u_a|/2) − G_ant|`, maximum over `p`:

  | `q` | `Δ` | `q²Δ` |
  |---|---|---|
  | 59 | `5.7e-4` | 1.98 |
  | 127 | `1.3e-4` | 2.10 |
  | 251 | `5.5e-5` | 3.47 |

  The `q²Δ` column is not constant, so these rows do not support the reading `O(q⁻²)` of the original. Only the table is claimed. The tail satisfies `|r_k| < 1e-16` at `k = 16`.
- **V15 (`κ` and `G` on four `q`), form T.** Over all `p` and `q ∈ {59, 127, 251, 1009}`: `Re κ ∈ [0.015, 0.064]`, `|Im κ| ≤ 0.055`, `|r₃| ≤ 2.3e-3`, `|r₄| ≤ 1e-4`. `κ(−x̃) = conj κ(x̃)` holds exactly by conjugation (form E). Newton on the cubic `1 − u + κu² + r₃u³ = −1` from `u = 2` reproduces `G` to `≤ 3e-3` on the same sets.
- **V16 (bounded `p`), form T.** Values at the largest computed `N`, and the last doubling step (`data/sequences.json`):

  | family | `p` | `N` | `q` | `G` | `κ` | `|G(N) − G(N/2)|` |
  |---|---|---|---|---|---|---|
  | `[0;3,N]` | 3 | 384 | 1153 | 1.12837 | `0.05453 + 0.00506i` | `1e-5` |
  | `[0;2,1,N]` | 3 | 384 | 1154 | 1.12537 | `0.05504 − 0.02351i` | `1e-5` |
  | `[0;2,N]` | 2 | 384 | 769 | 1.11142 | `0.04969 + 0.02399i` | `1e-5` |
  | `[0;N]` | 1 | 384 | 384 | 1.02632 | `0.02383 − 0.05270i` | `7e-5` |

  In the `[0;N]` row, `Im κ` still moves by `3.5e-4` between `N = 192` and `N = 384`.
- **V16b (generic tails), forms T + S.** `x̃ = [0; tail, N, a₁]`, `a₁ ∈ {16, 64}`, values at `N = 64, 128, 256` (`data/generic_limits.json`). `G₀` is the fitted intercept of `G = G₀ + c/N` on these three nodes (form S).

  | side | `G(64)` | `G(128)` | `G(256)` | `G₀` (S) |
  |---|---|---|---|---|
  | `⅓⁻` | 1.12345 | 1.12196 | 1.12122 | 1.1205 |
  | `⅓⁺` | 1.12152 | 1.12002 | 1.11927 | 1.1185 |
  | `½⁻` | 1.10505 | 1.10257 | 1.10133 | 1.100 |
  | `0⁺` | 1.01756 | 1.01293 | 1.01067 | `1.008 ± 0.002` |

  At `q = 16385` the `0⁺` value is `1.00998`. Between `a₁ = 16` and `a₁ = 64`, `G` differs by `< 3e-4`.
- **V17 (the far digit), form T.** At fixed `N`, bulbs whose `x̃` differ only in the deepest digit `a₁` (`data/limit_grid.json`):

  | `N` | `x̃ = [0;3,N]` (`p = 3`) | `x̃ = [0;3,N,2]` | `x̃ = [0;3,N,64]` | gap, `p = 3` vs `a₁ = 64` |
  |---|---|---|---|---|
  | 32 | 1.12792 | 1.12855 | 1.12614 | 0.00178 |
  | 64 | 1.12830 | 1.12688 | 1.12322 | 0.00508 |
  | 128 | 1.12835 | 1.12596 | 1.12173 | 0.00662 |

  At every tabulated `N`, `G` is monotone in `a₁` over `a₁ ∈ {2, 3, 5, 9, 16, 32, 64}`. The gap grows with `N` on all three rows.
- **V18 (arithmetic in `t`), form T.** The prefix table at tail `(N = 128, 3)` and the dense `q = 1009` decomposition stand as in the original. They are finite regressions on stated bins: 100 `x*`-bins leave `R² = 0.9984`, and adding 50 `t`-bins leaves `R² = 0.9989`.
- **V19 (component comparison), form T.** As in the original. The `disk2/main2` ratios are values at the stated `a₁` and `N`.
- **V20 (golden tail), forms T + S.** `κ = 0.05716, 0.06216, 0.06390` at `q = 34, 89, 610`. A geometric fit on these three nodes gives the statistic `κ₀ ≈ 0.0648 − 0.0058i`, `G₀ ≈ 1.158`.
- **V21 (`Im κ` steps), form T.** `J` is already finite: a difference of window means over `0.002 < |x̃ − p'/q'| < 0.009` at `q = 1009`. The original table of `J` and `J·q'²` stands unchanged.
- **V22 (P6 evaluated, `q ≤ 20`), forms E + T.** `bulbford/cycles.py` pulls back all `2^q − 1` rays of period dividing `q` at the root of `B_{p/q}`. The `q` rays of the `p/q` rotation cycle land at `z₀`, and the other `2^q − q − 1` land one each on the other fixed points of `F = f^q`, which are summed. Runs cover every `p ≤ q/2` with `2 ≤ q ≤ 20` (64 fractions; the rest are complex conjugates). Data: `data/cycle_index_sum.json`, `scripts/cycle_index_sweep.py`.
  (a) *T.* The cycle sum equals the Arb ball of `ι` to `≤ 1.8e-14` on all 64. For `q ≤ 12` the tests also check that the points are distinct and satisfy `f(z_j) = z_{2j}`.
  (b) *E.* Let `α_0 < … < α_{q−1}` be the rotation-cycle angles, `(α_w, α_{w+1})` the characteristic arc, and `M = 2^q − 1`. For every `p` and `3 ≤ q ≤ 20`, the doubling orbit of `α_{w−1} + 1/M` contains `α_{w+2} − 1/M` and has period `q`. This is the **flank cycle** (`flank_angles`). P7 proves it for every `q ≥ 3`: `w = p − 1`, and `α_{w+2} − 1/M` is `p̄` doublings after `α_{w−1} + 1/M`.
  (c) *T.* The flank cycle carries the largest `|term|` for all 63 fractions with `q ≥ 3`. Its share of `Σ|term|` is `0.04–0.79` (`0.04–0.11` at `q ≥ 17`), and its `|F′| ∈ [14.9, 51.3]`. For `q ≥ 4` it is not a rotation cycle; at `q = 3` it is the `2/3` rotation cycle.
  (d) *T.* For `p ≥ 2` at fixed `q`, the flank `|F′|` is strictly decreasing in `x* = |p̄/q|`: 30 of 30 adjacent comparisons for `q = 5 … 20`. At `q = 19` it runs from `51.29` (`x* = 2/19`) to `33.92` (`x* = 9/19`). `p = 1` is off this order (`47.93` at `q = 19`), the bounded-`p` exception of V16. For `p = 1`, the values for `q = 2 … 20` are `9.00, 14.93, 20.22, …, 47.93, 48.31`, with increasing values and decreasing increments.
  (e) *T.* The sum cancels heavily: `Σ|term| / |ι| ∈ [1.18, 7.40]`, and `5.5–7.4` for `p ≥ 2` at `q = 19, 20`. There the ten largest terms carry `0.12–0.20` of `Σ|term|`.
  (f) *T.* The period-`q` rotation cycles carry at most `0.084` of `Σ|term|` for `q ≥ 11` (largest at `p = 1`, `≤ 0.001` for `p ≥ 2` at `q = 20`). `β` carries `0.035–0.047` for `p = 1` at `q ≥ 8`, and `≤ 0.008` for `p ≥ 2`.
- **V23 (flank against `κ` and `G`, `q ≤ 20`), form T.** Over the 78 within-`q` pairs with `p ≥ 2` (`q = 5 … 20`), the flank `|F′|` is anti-ordered with `Re κ` on 67 pairs and with `G_ant` on 68. Every exception involves `p = 2`, the bulb with `x*` nearest `½`, where `Re κ` and `G` have passed their V15 peak near `x* ≈ 0.40` while the flank keeps falling (strictly decreasing in `x*`, V22(d)). At `q = 19`, `x* = 0.105 … 0.421` gives flank `51.29 … 34.91`, `Re κ` `0.0309 … 0.0490` and `G` `1.0560 … 1.1188`, while `x* = 0.474` gives flank `33.92`, `Re κ = 0.0443` and `G = 1.1031` (`scripts/flank_vs_kappa.py`, `data/flank_vs_kappa.json`).
- **V24 (the second ring, `q ≤ 20`), forms E + T (`scripts/second_ring.py`, `data/second_ring.json`).** Ring `k` around the characteristic arc is the pair `α_{w−k} + 1/M`, `α_{w+1+k} − 1/M`; ring 1 is the flank cycle (P7). (a) *Exact:* the two ring-2 angles lie on **distinct** doubling orbits, each of period `q`, for every `p` and every `q ≥ 5` (proven as P8). (b) *Numerical:* for `p ≤ q/2`, `5 ≤ q ≤ 20` (61 fractions), ranks 2 and 3 of the P6 sum are exactly these two cycles in 44. The 17 exceptions are exactly `p = 1` (`q ≥ 8`) and `p̄ = q − 2` (`q ≥ 13`), the two smallest values `x* = 1/q, 2/q`. There, a cycle through a neighbour of the arc itself intrudes at rank 2 or 3: `θ₊ + 1/M` for `p̄ = q − 2`, and `α_{w−1} − 1/M` or `θ₋ − 1/M` for `p = 1`. The ring-2 cycles then sit at ranks 2–5.
- **V25 (the small-`x*` intruders are the widest arcs), form T (`scripts/second_ring.py`, `data/second_ring.json`).** For `p ≤ q/2` and `5 ≤ q ≤ 20`, all 23 cycles at rank 2 or 3 that are not ring-2 cycles pass through an endpoint of one of the two widest arcs (widths `2^{q−1}/M`, `2^{q−2}/M`). The V24 exceptions all have `x* ≤ 2/q`, that is, a flanking arc is among the two widest (P9). The members of that family that are not exceptions are exactly those below the thresholds of V24(b): `p = 1` with `q ≤ 7`, and `p = (q−1)/2` with `q ≤ 11`. **Reading:** V24's "small-`x*` effect" is a geometric fact about rays. `x*` is the exponent, over `q`, of the narrower arc next to the characteristic arc. When `x* ≤ 2/q`, the other flanking arc is a quarter or half of the circle, and a cycle through that arc's endpoint competes with ring 2. What is still numerical is only the crossing, at `q = 8` for `p = 1` and `q = 13` for `p = (q−1)/2`.
- **V26 (the crossings), form T (`scripts/crossings.py`, `data/crossings.json`).** For the two families of V25, `R(q) = max |intruder term| / min |ring-2 term|` is computed from single orbits (`cycle_through`, verified by the ray relation and `|F(z) − z|`) for `q ≤ 64`. Here the intruders are the cycles through widest-arc endpoints other than rings 1 and 2. `R > 1` agrees exactly with the V24 exceptions for `q ≤ 20`. **`p = 1`:** `R(7) = 0.943`, `R(8) = 1.010`. `R` increases at every step for `q = 6 … 64`, reaching `4.239`. **`p = (q−1)/2`:** `R(11) = 0.932`, `R(13) = 1.012`. `R` increases at every step for `q = 9 … 63`, reaching `2.437`. Both crossings are narrow, about `1%`. All terms grow with `q`, and the leading intruder grows fastest: for `p = 1` it goes `0.0801 … 1.1831` over `q = 8 … 64`, against `0.0794 … 0.2791` for the smaller ring-2 term and `0.2310 … 1.2585` for the flank.
- **V27 (term growth), forms T + S (`scripts/multiplier_growth.py`, `data/multiplier_growth.json`).** A period-`q` cycle contributes `q/|1 − F'|` in size, so its growth in `q` is read from `|1 − F'|`. Single verified orbits give `|1 − F'|` for `q = 8, 16, …, 256` (`p = 1`) and `q = 9, 17, …, 249` (`p = (q−1)/2`). The fit `a + b/q + c/q²` on the nodes `q ≥ 64` gives these intercepts `a` (form S):

  | cycle | `p = 1` | `p = (q−1)/2` |
  |---|---|---|
  | flank | `51.18` | `56.21` |
  | ring 2, lower / upper | `175.78` / `233.42` | `260.50` / `228.54` |
  | widest-arc intruder (`α_{w−1} − 1/M` for `p = 1`; `θ₊ + 1/M` for `(q−1)/2`) | `80.62` | `84.00` |
  | `β`, as `q·\|1 − F'_β\|` (P11) | `39.499` | not bounded: `\|2 − λ₀\|` is near `3` |

  The largest residual is `6.2e-03` for `p = 1` and `2.3e-02` for `p = (q−1)/2`. Every named cycle except `β` has a multiplier that settles, so its term grows linearly in `q`. For `p = 1` the leading intruder is `α_{w−1} − 1/M` for `q ≤ 21` and `β` for every `22 ≤ q ≤ 64`, because `β`'s slope `1/39.50` exceeds `1/80.62`. V26's `R` uses the smaller ring-2 term, which is the one with the larger `|1 − F'|`: ring-2 upper (`233.42`) for `p = 1` and ring-2 lower (`260.50`) for `p = (q−1)/2`. The ratios of the fitted intercepts, `233.42/39.50 = 5.91` (`p = 1`) and `260.50/84.00 = 3.10` (`p = (q−1)/2`), are the fitted plateaus of V26's `R(q)`. The largest computed values are `R(64) = 4.239` and `R(63) = 2.437`.

---

## 4. FALSIFIED

- F1–F3: as in the original.
- **F4 (no Lipschitz function of `x̃`), form X.** Any function `g` of `x̃` with `G = g(x̃)` on the pairs below must have a difference quotient `|ΔG/Δx̃|` at least as large as the last column. `x̃` is exact.

  | bulb A (`G`) | bulb B (`G`) | `|Δx̃|` | `|ΔG|` | slope |
  |---|---|---|---|---|
  | `1/256` (1.02628) | `16129/16385` (1.00998) | `2.38e-7` | 0.01631 | `6.8e4` |
  | `766/769` (1.12836) | `769/12307` (1.12122) | `1.06e-7` | 0.00714 | `6.8e4` |
  | `382/385` (1.12835) | `385/24643` (1.12173) | `1.05e-7` | 0.00662 | `6.3e4` |
  | `3/770` (1.12537) | `11553/12323` (1.11927) | `1.05e-7` | 0.00610 | `5.8e4` |
  | `511/513` (1.11141) | `513/8210` (1.10133) | `2.37e-7` | 0.01008 | `4.2e4` |

  Across all 504 bulbs at `q = 1009`, `G` ranges over `[1.0144, 1.1571]` (width `0.143`). The `x*`-bin fit of V18 leaves residual sd `1.6e-3`. So `G = g(x̃)` with `g` of Lipschitz constant `< 6.8·10⁴` is false on this finite data. Every `|ΔG|` exceeds the `1e-5` doubling stability of V16 by a factor of at least 600. This is the finite content of the original F4. The five pairs have the largest slopes among the pairs with `|ΔG| ≥ 0.005`, found by scanning all 678 computed bulbs.
- **F5 (the `(q_{n−1}/q)²` weight), form X.** The model "the digit `a₁` enters with weight `≤ C/N²`" requires `C ≥ gap · N²` on the rows of V17: `C ≥ 1.8, 20.8, 108` at `N = 32, 64, 128`. The required constant grows with `N` on the computed rows, and the model with any `C < 108` is false at `N = 128`.
- **F6 (the Farey-neighbour prediction of the P6 move), form X.** Prediction: "the cycles of the Farey neighbours `|pq' − p'q| = 1` carry the largest share". For prime `q`, a neighbour with `1 < q' < q` has `q' ∤ q`, so its cycle is not a fixed point of `f^q` and has no term in the sum. The one neighbour with a term, `β` (`q' = 1`), carries `≤ 0.047` for `q ≥ 6`, and the period-`q` rotation cycles carry `≤ 0.084` (V22(f)). The largest term belongs to the flank cycle in all 63 fractions (V22(c)).
- ~~V6's equal one-sided values at `⅓`~~ superseded, as in the original.

---

## 5. CLASSICAL REFERENTS (never asserted here)

Each row is a statement about infinitely many bulbs. It is recorded so that the finite work has a named target, together with the finite shadow that every new computation can check.

| id | classical statement (not asserted) | finite shadow checked so far | status |
|---|---|---|---|
| [GM84] | `diam B_{p/q} = O(q⁻²)` | `G ∈ [1.0144, 1.1571]` on all computed `p` at `q = 1009`; `[1.0231, 1.1559]` at `q = 251` | imported, unused |
| C1‴ | `G_H = Ĝ(x̃) + Ĥ_H(p; x̃) + o(1)`, `Ĝ` universal | V18 regression (`R² = 0.9984` from `x*` alone); V19 ratio `1.0000 ± 0.0002` for `a₁ ≥ 2` | conjectured |
| C14″ | `κ̂(x̃)` exists, cusps in `Re`, jumps in `Im` at rationals | V21 steps `J` at `q' ≤ 7`; `J·q'² ∈ [−0.167, −0.106]` for `q' ≤ 6` | conjectured |
| C16′ | two-ended resonance sum for `ι/q` | V16 vs V16b gaps at `p = 1, 2, 3` | conjectured |
| C17′ | `ι_{p/q} = ½ + qκ̂(x̃) + O(1)` | P6 gives `ι` as a finite exact sum; the dominant terms are untested | conjectured |
| P11′ | `q·\|1 − F'_β\|` tends to `4π²` for `p = 1` (proven in the original register) | exact values `54.096` (`q = 64`), `42.668` (`q = 256`), `40.248` (`q = 1024`), `39.669` (`q = 4096`) |  proven, classical |
| C20 | `κ(1/q)` has a Lavaurs-phase value as `q` grows | V16 row `[0;N]` at `N ≤ 384` | untested |
| C2′, C9′, C8, C11, C15′ | as in the original | none new | unchanged |

---

## 6. Relationship to Ford circles (finite statement)

| level | Ford circles | bulbs |
|---|---|---|
| combinatorics | Farey adjacency and mediants, exact | identical and exact for `q ≤ 16` (P4) |
| leading scale | radius `1/(2q²)`, exact | `q²·d/(2|c_H'|) = G ∈ [1.0144, 1.1571]` on every computed bulb at `q = 1009` (T); [GM84] only in §5 |
| first correction | none: the `PSL(2,ℤ)` action on `ℙ¹(ℚ)` is exact | `[ε²]ρ = q³(ι − ½)` with `ι ∈ ℚ(ζ_q)` (P2, E); `ι` exact for `q ≤ 8` (P5) |
| symmetry | `PSL(2,ℤ)` on horoballs, exact | Galois `ζ_q ↦ ζ_q^{p'p̄}` on `ι` (E); `p ↦ p̄` organises V15–V21 (T) |

---

## 7. Next moves (finite forms)

1. **The flank cycle.** P6 is evaluated for `q ≤ 20` (V22), its Farey prediction is falsified (F6), and the flank term is compared with `κ` and `G` (V23). V22(b) is proven as P7, and the second-rank terms are the second ring except at the two smallest `x*` (V24). V24(a) is P8, and the small-`x*` intruders are the widest arcs next to the characteristic arc (P9, V25). The crossings are located and `R` increases after them to `q = 64` (V26). Term growth: `β` in closed form (P11), and the other named cycles have settling multipliers (V27). Next: relate the multiplier plateaus to C20 (form T).
2. **Steps and cusps.** Tabulate `J(p'/q')` and the cusp depths at `q = 2003` for `q' ≤ 7`, next to the `q = 1009` table (form T).
3. ~~**Certified antipodes.**~~ Done as P10 (`q ≤ 16`, form C). Next: a certified replacement for `SatelliteLabel`, and boxes at `q = 59` beside the V14 table.
4. **`z³ + c`.** Repeat V15–V19 as tables (`MAIN3` is wired).

---

## 8. Map from the original register

| original | here | change |
|---|---|---|
| P0, P2 | P0 + P2 | one Laurent identity; "total index converges" becomes `[ε⁰]I = ι` |
| P1, P3, P4, P5 | P1, P3, P4, P5 | already finite; unchanged |
| none | P6 | projective index sum, exact at `ε = 0` |
| P10 | P10 | certified antipodes and `G_ant` brackets; unchanged (form C) |
| V12, V13, V15, V19, V21 | same ids | remainders and arrows removed; values unchanged |
| V14 | V14 | table only; the `O(q⁻²)` reading is withdrawn (`q²Δ` = 1.98, 2.10, 3.47) |
| V16, V17, V20 | V16, V16b, V17, V20 | limits become values at the largest `N` (T) or fit intercepts (S) |
| V18 | V18 | regressions on stated bins |
| F4 | F4 | "no single continuous `Ĝ`" becomes a Lipschitz lower bound `6.8·10⁴` from exact pairs |
| F5 | F5 | "survives `N → ∞`" becomes `C ≥ 1.8, 20.8, 108` at `N = 32, 64, 128` |
| C1‴, C14″, C16′, C17′, C20, C2′, C9′, C8, C11, C15′ | §5 | classical referents with finite shadows |
| §6 "asymptotic at leading order" | §6 | a two-sided bound on computed bulbs |
