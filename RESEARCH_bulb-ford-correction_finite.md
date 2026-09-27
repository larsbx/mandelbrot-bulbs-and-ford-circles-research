# RESEARCH BRANCH, finite version: bulb geometry of M vs Ford circles, without limits

Branch id: `bulb-ford-correction`  ·  finite revision 2026-09-27 of `RESEARCH_bulb-ford-correction.md` (2026-09-25)  ·  status: **active**
Discipline: PROVEN / VALIDATED / CONJECTURED / FALSIFIED are absolute, as in the original register. This version adds a second rule. **No claim below is a limit.** Every limit, `O(·)`, `→` and "converges" in the original is replaced by one of the statement forms in §0, or moved to §5 as a named classical referent that this document never asserts. The original register is unchanged and remains the claim-state file of record. §8 maps every original item to its form here.

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
  as an exact identity at `ε = 0`, over the finitely many other fixed points. This is the limit-free form of the route in C17′: it evaluates `ι` as a finite sum rather than as the index of a merger.

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

  The `q²Δ` column is not constant, so these three rows do not support the original's reading "i.e. `O(q⁻²)`". Only the table is claimed. The tail satisfies `|r_k| < 1e-16` at `k = 16`.
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

1. **P6 at work.** Evaluate `ι_{p/q} = −Σ 1/(1−ρ_z)` exactly for `q ≤ 20` over the dynatomic factors, and tabulate the share of each cycle. The expectation to test is that the cycles of the Farey neighbours `|pq' − p'q| = 1` carry the largest share (form T).
2. **Steps and cusps.** Tabulate `J(p'/q')` and the cusp depths at `q = 2003` for `q' ≤ 7`, next to the `q = 1009` table (form T).
3. **Certified antipodes.** A two-variable Krawczyk witness for `(f_c^q(z) − z, (f_c^q)'(z) + 1)`, turning `G_ant` for `q ≤ 16` into form C.
4. **`z³ + c`.** Repeat V15–V19 as tables (`MAIN3` is wired).
5. **An audit.** A prose check that rejects `→`, `lim`, `O(` and "converges" outside §0, §5 and §8, in the manner of the terminology audits of the sibling repositories.

---

## 8. Map from the original register

| original | here | change |
|---|---|---|
| P0, P2 | P0 + P2 | one Laurent identity; "total index converges" becomes `[ε⁰]I = ι` |
| P1, P3, P4, P5 | P1, P3, P4, P5 | already finite; unchanged |
| none | P6 | projective index sum, exact at `ε = 0` |
| V12, V13, V15, V19, V21 | same ids | remainders and arrows removed; values unchanged |
| V14 | V14 | table only; the `O(q⁻²)` reading is withdrawn (`q²Δ` = 1.98, 2.10, 3.47) |
| V16, V17, V20 | V16, V16b, V17, V20 | limits become values at the largest `N` (T) or fit intercepts (S) |
| V18 | V18 | regressions on stated bins |
| F4 | F4 | "no single continuous `Ĝ`" becomes a Lipschitz lower bound `6.8·10⁴` from exact pairs |
| F5 | F5 | "survives `N → ∞`" becomes `C ≥ 1.8, 20.8, 108` at `N = 32, 64, 128` |
| C1‴, C14″, C16′, C17′, C20, C2′, C9′, C8, C11, C15′ | §5 | classical referents with finite shadows |
| §6 "asymptotic at leading order" | §6 | a two-sided bound on computed bulbs |
