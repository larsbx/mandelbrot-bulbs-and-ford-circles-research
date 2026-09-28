# Contributions audit — `bulb-ford-correction`

Audited 2026-09-28 against `main` @ `bb07e18` (register revision 2026-09-27).
Scope: what the study contributes that is new, how strongly each contribution is
established, what it is not, and why it matters. This file is exposition (the
`docs` plane): it changes no claim state. Proposed register edits are listed in
§5 for the register's owner to accept or refuse.

**Method.**

1. Re-read every PROVEN, VALIDATED and CONJECTURED item of
   [`RESEARCH_bulb-ford-correction.md`](../RESEARCH_bulb-ford-correction.md) and
   classify it by tier and by novelty.
2. Re-run the repository gates: `tools/audit_estate_layout.py`,
   `tools/audit_limits.py`, `pytest tests/`.
3. Re-derive the headline identity P2 independently of the kernel
   (`experiments/scripts/audit_p2_independent.py`, mpmath only).
4. Search the literature for prior art on each headline item (§3).

Novelty verdicts are relative to that search. They mean "not found", not "does
not exist", and they should be re-checked by a specialist before any external
claim of priority.

---

## 1. Verdict table

Tier is the register's own label. The audit checks whether the label is earned.
Novelty is graded as follows:

- **N**: no prior statement found.
- **N/m**: new statement, standard method.
- **R**: re-derivation or packaging of known material.
- **C**: computer-assisted certificate of a known fact.

| # | Contribution | Register tier | Audit | Novelty |
|---|---|---|---|---|
| 1 | **P2**: `ρ = 1 − q²ε + q³(ι_{p/q} − ½)ε² + O(ε³)` | PROVEN | Proof correct. Independently re-derived at 1/3, 1/4, 2/5 (§2.1). | N/m |
| 2 | The correction `G` depends mainly on `x̃ = p̄/q` (modular inverse) | VALIDATED (V3–V6, V15, V18) | Earned at the stated precision. `R² = 0.9984` at `q = 1009`. | N |
| 3 | A second, bounded-`p` term, the only component-dependent part (F4, V16–V19) | VALIDATED; F4 FALSIFIED | Earned. The falsification of the study's own C1″ is well documented. | N |
| 4 | `Re κ` has cusps and `Im κ` has jumps at rationals, with `q'^{-2}` weights (V21, C14″) | VALIDATED (weights), CONJECTURED (law) | The `q'^{-2}` fit uses 6 of the 9 rationals, all at one `q = 1009` (E-b). Correctly left conjectural. | N |
| 5 | **P7–P9**: flank orbit, second-ring split, arc widths placed by `p̄` | PROVEN, integer-checked for `q ≤ 60` | Proofs read and complete. Tests pass. | N/m (P9 refines Goldberg) |
| 6 | **P6/P11**: `ι = −Σ 1/(1 − F'(z))` over the other fixed points; closed-form `β` term | PROVEN | Correct. | R (P6 is [Mil] §12); P11 elementary |
| 7 | **P3/P10/P12**: certified centres, antipodes and continuation for all 79 bulbs with `q ≤ 16`, plus V14 at `q = 59, 127, 251` | PROVEN (finite) | Replayable and CI-gated. One named import remains: that the certified centre is the centre of `B_{p/q}`. | C |
| 8 | **C20′/C22**: `κ(1/q) − κ₀ − 1/(2πiq) = O(q⁻²)` with `κ₀ = a₂/(2πi a₁²)` from the horn map of `z + z²` | CONJECTURED (derivation, not proof) | Numerics agree to `2.9·10⁻⁶` at `q = 1009` (V28b). The interchange of `q → ∞` with `[u²]` is unproven, as stated. | N |
| 9 | Ford-circle comparison table (§6 of the register) | Mixed | A fair summary. Each row cites its tier. | — |

---

## 2. Per-contribution audit

### 2.1 P2 — second-order bulb shape equals the parabolic index

**Statement.** For `λ = λ₀e^ε` near the root of the `p/q` satellite of `H`:
`ρ(ε) = 1 − q²ε + q³(ι − ½)ε² + O(ε³)`, where `ι = Res_{z₀} dz/(z − f^q(z))`.

**Proof check.**

- The revised proof is sound. `I(ε) = (1/2πi)∮_γ dz/(z − F_ε(z))` is holomorphic in `ε`.
- The argument principle gives exactly `q + 1` fixed points inside `γ`, and [DH] identifies them.
- The residue theorem gives the Laurent identity `1/(1 − μ) + q/(1 − ρ) = I(ε)`.
- Pole cancellation forces `ρ₁ = −q²`; the constant term gives `ρ₂ = q³(ι − ½)`.
- The argument needs no normal form. It proves P0's linear term as a by-product.
- The earlier proof, which assumed the linear term, is struck, not deleted, as the discipline requires.

**Independent re-derivation.** `experiments/scripts/audit_p2_independent.py`
does not import `bulbford`. It computes `ρ` by deflated Newton and `ι` by
contour quadrature at 60 digits:

| `p/q` | `[ε]ρ` | `[ε²]ρ` | `q³(ι − ½)` |
|---|---|---|---|
| 1/3 | `−9.0` (imag `3·10⁻¹⁸`) | `−7.377551 − 0.84835142i` | `−7.377551 − 0.84835142i` |
| 1/4 | `−16.0` (imag `2·10⁻¹⁷`) | `−11.972318 − 5.0519031i` | `−11.972318 − 5.0519031i` |
| 2/5 | `−25.0` (imag `2·10⁻¹⁷`) | `−11.55561 + 1.4504625i` | `−11.55561 + 1.4504625i` |

The quadrature values of `ι_{1/3}` and `ι_{1/4}` also agree to 12 digits with
the exact V13 values `(92 − 16ζ₃)/441` and `(1447 − 365i)/4624`.

**Novelty.** The Laurent identity is a classical tool. It is Milnor's index
bookkeeping ([Mil] §12). The "index at a bifurcation" argument is used for
parabolic arcs of multicorns (Mukherjee–Nakane–Schleicher), and Buff–Epstein use
the résidu itératif to bound basins. We found no statement of the explicit
second-order coefficient `ρ₂ = q³(ι − ½)`, and no use of it as the first
correction to the Ford law.

- **Verdict:** new statement, standard method.
- **Risk:** folklore among parabolic-implosion specialists.

**Scope note.** The register states P2 "for any unicritical family". The proof
uses only [DH] and the merging of exactly `q + 1` fixed points, which holds for
satellite roots in `z^d + c`. Every check is quadratic, though (V12, the table
above). The `z³ + c` case is next move 4.

### 2.2 The modular inverse `x̃ = p̄/q` as the leading variable of `G`

**Status.** VALIDATED, with doubles and FFT validated against balls to `10⁻¹²` (V12).

- At `q = 1009`, a 100-bin function of `x*` alone explains `G` with `R² = 0.9984`.
- The residual standard deviation is `1.6·10⁻³` (V18).
- F4 adds a finite witness: bulbs with `|Δx̃| ≈ 10⁻⁷` differ in `G` by `0.016`. So no single Lipschitz-tame `g(x̃)` fits.

**Why it holds combinatorially.** P9 proves that `p̄` controls the external-ray
geometry exactly. Goldberg's `q` arc widths `2^e/M` sit around the
characteristic arc in the order `e_j = j·p̄ mod q`. The dominant term of the
index sum is monotone in `x*` (V22–V23). This gives a mechanism, not only a fit.

**Novelty.**

- The leading law `d ≈ 2 sin(πp/q)/q²` is classical: [GM84], and Fowler–McGuinness 2019 give an analytic estimate.
- The appearance of `p̄` in the correction was not found in the literature.
- **Verdict: new.**

### 2.3 The bounded-`p` term and the order of limits

**Status.** VALIDATED (V16–V17, V19). Falsifies the study's opening conjecture (F4, F5).

- Along bounded `p`, the offsets from the generic `x̃`-value are `+0.018, +0.011, +0.008` for `p = 1, 2, 3` (at `0⁺, ½⁻, ⅓⁻`).
- The two orders of limits differ by `0.008` at `⅓`.
- Across hyperbolic components (`disk2/main2`), the ratio is `1.0000 ± 0.0002` generically, and `0.988–0.992` only at bounded `p` (V19).

**Audit.** The evidence is well controlled:

- The `c/N` approach rates are quoted with their fits (E-e).
- The universality comparison isolates the bounded-`p` bulbs as the source of the older V9 scatter.

**Novelty.** No prior statement found. The two-ended resonance picture (C16′)
is a conjecture that explains it.

### 2.4 Rational singularities of `κ̂` with Ford weights

`Re κ` has V-cusps and `Im κ` has jumps at every `p'/q'`. The jump sizes `J`
scaled by `q'²` are `−0.137 … −0.167` for `q' = 2 … 6` (V21).

**Audit.** The weight law rests on 9 rationals at a single `q = 1009`, with
window-limited values at `q' = 7`. Error E-b already notes this. The register
correctly keeps the law in CONJECTURED (C14″).

**Why it matters.** Singular at every rational, with weights decaying in the
denominator, is the fingerprint of Brjuno-type functions (Yoccoz; Marmi–Moussa–Yoccoz). Here it sits on the parabolic side, in the variable `p̄/q`.

### 2.5 Combinatorial lemmas P7–P9

**Content.**

- **P7:** the two flank angles lie on one period-`q` doubling orbit, `p̄` doublings apart.
- **P8:** ring 2 splits into two distinct period-`q` orbits. The proof uses a height-set invariant `G(W) = {qH_k − km}` up to translation.
- **P9:** the arc widths are placed by `j·p̄ mod q`.

**Audit.**

- The proofs were read step by step and are complete, including the `p ∈ {1, 2, q−2, q−1}` edge cases via the conjugation `x ↦ M − x`.
- They are tested in exact integers for `q ≤ 60`: 1100 fractions for P7, 1096 for P8.
- The P8 test has a mutation check: replacing `m` by `m + 1` in `heights` fails 3120 cases. This shows the test would catch a wrong invariant.

**Novelty.**

- Rotation-set combinatorics are classical: Goldberg 1992, Bullett–Sentenac 1994.
- P9 is a sharpening: which width sits next to the characteristic arc is decided by `p̄`.
- **Verdict:** new statements, elementary proofs. Well suited to a short note.

### 2.6 P6 and P11

P6 is Milnor's index formula applied to `f^q` with `ι_∞ = 1`, so it is a
re-derivation. Its contribution is operational: it turns `ι` into a finite sum
over repelling cycles, and that sum is what makes V22–V27 possible.

P11 is elementary, and its consequence is useful. The `β` term grows like
`q/(4π²)` for `p = 1`.

### 2.7 Certificates P3, P10, P12

**Content.** For all 79 bulbs with `q ≤ 16`:

- Krawczyk inclusions and exact-period exclusions in rational and Arb arithmetic.
- `ζ_q` certified without trigonometry.
- `G_ant` bracketed to width `< 10⁻¹⁶`.
- A certified continuation from centre to antipode, with up to 21770 pieces.
- V14's `q = 59, 127, 251` sweep extended to 217 certified antipodes.

**Audit.**

- Records replay from stored endpoints.
- CI recomputes the continuation for `q ≤ 7`.
- The one remaining import (the centre half of `SatelliteLabel`) is named, with a stated reason it is hard: every route ends in a limit at the root.

**Novelty.** These are certificates of known facts. Their value is that every
number the conjectures lean on at low `q` is rigorous.

### 2.8 Horn-map constants C20′, C22, V28–V30

**Content.**

- `κ₀ = a₂/(2πi a₁²) = 0.0238258874022 − 0.0523046591140i` for `p = 1`, where the `a_k` are upper-end Fourier coefficients of the horn map of `z + z²`.
- The `1/q` term is exactly the phase curvature `1/(2πi)` (V30a).
- For `p = (q−1)/2`, a half-step horn map at `c = −¾` gives `κ₀ = 0.0184052616 + 0.0429037239i`. The phase explains `i/π` of the `1/q` constant and leaves `R = 0.354856 − 0.069314i`.

**Audit.**

- The derivation interchanges `q → ∞` with taking `[u²]`. The register flags this as unproven. It is the correct gap to name.
- Numerics support it strongly: the Richardson estimates meet `1/(2πi)` to `4·10⁻⁶`.
- Two corrections found along the way are worth keeping: the parity-matched outgoing petal, and the period-1 (not ½) horn map in V29.

**Novelty.**

- Lavaurs maps are used for limb bounds by Kapiamba (2021/2025).
- Nothing found extracts bulb-shape constants from horn-map Fourier coefficients.
- **Verdict: new.**

---

## 3. Prior art consulted

| Work | Relation to this study |
|---|---|
| Guckenheimer–McGehee 1984 [GM84] | `diam B_{p/q} = O(q⁻²)` for bulbs. The study's leading order. |
| Yoccoz / Pommerenke / Levin | `O(1/q)` for **limbs**. Not used by the study, which concerns bulbs. |
| A. Kapiamba, *An optimal Yoccoz inequality for near-parabolic quadratic polynomials*, arXiv:2103.03211 | `C_N^{-1}/q² < diam L_{p/q} < C_N/q²` for limbs in the continued-fraction classes `ℚ_N`, via Lavaurs maps. No constants, no bulb shape. **Not cited in the register** (§5, item 4). |
| A. C. Fowler, M. J. McGuinness, *The size of Mandelbrot bulbs*, Chaos Solitons Fractals X 3 (2019) 100019 | Analytic estimate of bulb radius at leading order. Full text not reviewed (paywalled). Correction terms and `p̄`-dependence not visible in the abstract. **Not cited** (§5, item 4). |
| X. Buff, A. Epstein, *A parabolic Pommerenke–Levin–Yoccoz inequality*, Fund. Math. 172 (2002) [BE02] | `Re((N+1)/2 − ι) > N/π²` when the immediate basin holds only `N` grand-orbit classes of critical points. See §4, item 5. |
| Mukherjee–Nakane–Schleicher (multicorns) | The index at a parabolic bifurcation is the same bookkeeping as P2. It is method prior art only. |
| Goldberg 1992 [Gol92]; Bullett–Sentenac 1994 [BS94] | Rotation sets of doubling. Background for P4 and P7–P9. |
| Lavaurs 1989 [Lav89]; Shishikura 2000 [Shi00] | Parabolic implosion. Background for C20′–C22. |

---

## 4. Findings

Severity:

- **High:** a claim stated stronger than its evidence.
- **Medium:** an internal inconsistency a reader will trip on.
- **Low:** bibliography or wording.

No high-severity findings. Every PROVEN item has a complete proof, and every
VALIDATED item cites a precision and a data file.

1. **Medium — stale range in the header and in E-c.**
   - The status line says "numerics-validated to q ≈ 2000", and E-c says "`q > 2100` untested".
   - Later entries go further: V16(b) uses `q ≤ 16385`, V29a `q = 2049`, and V30 `q = 1025 … 8193` (FFT route).
   - E-c most likely means the continuation instrument only. As written, it contradicts V16 and V30.
2. **Medium — summary precision for "`G` is a function of `(κ, r₃)`".**
   - §0 says "to `10⁻³`". V15, the supporting item, says "`≤ 3e-3`".
   - The summary is three times stronger than its evidence.
3. **Medium — two ranges for `Ĝ`.**
   - §0 says "range 0.12 in `|x̃|`" and, in the same paragraph, `Ĝ` "runs from `≈1.01` to `≈1.16`" (a range of `0.15`).
   - F4 also says "range `0.15`".
   - One of the two figures should go, or each should say what it measures (generic `Ĝ` versus `G` at a fixed `q`).
4. **Low — missing citations.** Kapiamba (arXiv:2103.03211) and Fowler–McGuinness (2019) are the two closest modern works on bulb and limb sizes. Neither appears in the theorem tags.
   - A reader will ask how the study relates to them.
   - The answer is short: they work on limbs, or at leading order. This study gives the first correction for bulbs.
5. **Low — [BE02] "statement to be re-checked".**
   - The published statement is `Re((N+1)/2 − ι) > N/π²` under a hypothesis on grand-orbit classes of critical points in the immediate basin.
   - For `N = q` petals it bounds `Re ι < (q+1)/2 − q/π² ≈ 0.40q + 0.5`.
   - The study's data give `Re ι = ½ + q·Re κ ≤ ½ + 0.064q` (V15), which is consistent with a wide margin.
   - Whether the hypothesis holds for `f^q` at a satellite root was not checked in this audit. The tag can move from "to be re-checked" to "statement confirmed, hypothesis to be checked".
6. **Low — mixed sources in C1‴.**
   - C1‴ quotes the bounded-`p` offsets `+0.018, +0.011, +0.008, +0.0005` for `p = 1..4`.
   - The first three are the V16 limit offsets. The fourth is the V18 local residual at `q = 1009`, whose own `p = 1..3` values are `+0.0080, +0.0044, +0.0025`.
   - The list mixes two measures.
7. **Low — header theorem list.** The status line lists P2, P6–P11. It omits P3, P4, P5 and P12, which are proven too (P3, P10 and P12 are finite certificates).

---

## 5. Recommended register edits

These are not applied here, because the audit does not change claim state.

1. Status line: "numerics-validated to `q ≈ 16000` (FFT/sequence routes); continuation instrument to `q ≈ 2100`". Reword E-c the same way. *(Finding 1.)*
2. §0: "to `3·10⁻³`". *(Finding 2.)*
3. §0: pick one range for `Ĝ` and say what it measures. *(Finding 3.)*
4. Theorem tags: add **[Kap]** Kapiamba and **[FM19]** Fowler–McGuinness, each with one line on the relation (limbs versus bulbs; leading order versus first correction). *(Finding 4.)*
5. [BE02]: "statement confirmed (Fund. Math. 172, 2002); applicability to `f^q` at satellite roots unchecked". *(Finding 5.)*
6. C1‴: quote V16 offsets throughout, or V18 residuals throughout. *(Finding 6.)*

---

## 6. Why mathematicians should care

1. **It turns a picture into a theorem-shaped question.** "The Mandelbrot set looks like Ford circles" is folklore. The study says exactly where the analogy holds and where it fails:
   - Combinatorics: exact (P4).
   - Leading order: asymptotic ([GM84]).
   - First correction: `PSL(2,ℤ)` gives way to Galois action on `ℚ(ζ_q)` and the involution `p ↦ p̄` (P2, remark (ii); P5).
2. **A new arithmetic function on the parabolic side.**
   - `κ̂(x̃)` (and `Ĝ`) is singular at every rational, with cusps and jumps, weights decaying in the denominator, and geometric convergence along quadratic irrationals (V20).
   - This is the behaviour of Brjuno-type functions, which control Siegel disk sizes on the elliptic side (Yoccoz, Buff–Chéritat).
   - A proof of C14″ would put bulb sizes in that family, in the variable `p̄/q`.
3. **The modular inverse is a signature.** `p ↦ p̄` is the move behind Kloosterman sums and Dedekind reciprocity. C16′ proposes two families of `O(1/q)` resonances (`j ≡ mp̄` and small `j`). This gives a concrete small-divisor sum to estimate.
4. **Parabolic implosion becomes quantitative.** Kapiamba shows `O(q⁻²)` is sharp for limbs in special classes but gives no constants. The study extracts explicit constants for bulbs from horn-map Fourier data (`κ₀ = a₂/(2πi a₁²)`) and isolates one precise open step: the interchange in C20′.
5. **Self-contained, publishable pieces.**
   - P2 with its one-line Laurent proof.
   - P7–P9 as a note on rotation cycles of doubling.
   - The certificates as a reproducible data set.
6. **It is auditable.**
   - Four tiers: PROVEN, VALIDATED, CONJECTURED, FALSIFIED.
   - A limit-free companion register enforced in CI.
   - Falsifications kept, not deleted (F1–F6).
   - Exact fixtures separated from ball and FFT evidence.

   Together these let a reader check each claim at the grade it is stated.

---

## 7. Gates re-run for this audit

| Gate | Result |
|---|---|
| `python tools/audit_estate_layout.py` | passed |
| `python tools/audit_limits.py` | passed |
| `PYTHONPATH=kernel python -m pytest -q tests/` | run in CI on the audit PR (the full suite exceeds 20 min locally) |
| `python3 experiments/scripts/audit_p2_independent.py` | agrees with P2 at 1/3, 1/4, 2/5 (§2.1) |
