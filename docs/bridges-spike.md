# Bridges spike — from bulb shapes to arithmetic

2026-09-28, on `main` @ `2b1b41a`. Exploratory: this file is exposition (the `docs` plane) and
changes no claim state itself; the register entries it supports were added after review (§6).

- **Script:** `experiments/scripts/bridges_spike.py` (2.6 s)
- **Data:** `experiments/data/bridges_spike.json`
- **Kernel:** `kernel/bulbford/norms.py`
- **Tests:** `tests/test_parabolic_norms.py` (57 tests)

The [contributions audit](contributions-audit.md) §6 named the bridges that would make the study
matter outside holomorphic dynamics. This spike tests four of them against data rather than
argument:

| # | Bridge | Question | Result | Tier |
|---|---|---|---|---|
| B1 | modular inverse | Is `p̄/q` really the variable, spectrally? | Yes. Fourier mass in `p̄/q` is 100× that in `p/q` at `m = 1`, and 10× in `ℓ²` over `m ≤ 6`. | VALIDATED (`q = 1009`) |
| B2 | Ramanujan sums / Ford weights | Does the spectrum carry the arithmetic of the rational singularities? | Yes for `Im κ`: a Ramanujan-sum expansion recovers the jump law at every denominator at once. No for `Re κ` under a V-cusp model. | VALIDATED, one `q` |
| B3 | wake combinatorics ↔ ℚ(ζ_q) | Is there arithmetic in `ι`'s denominators? | Yes: the primitive primes of `2^q − 1` and `4^q − 1` divide `N(a_q)`, with multiplicities predicted exactly. | PROVEN for `λ ∈ {2, 4}`; VALIDATED for `λ = −2` (27/27, `q ≤ 20`) |
| B4 | Dedekind sums | Is `s(p, q)` the two-ended variable of C16′? | No. It removes 3.9% of the residual variance (`1/p`: 1.2%), and most of that comes from the `p ≤ 4` bulbs. | negative |

---

## B1 — the spectrum lives in `p̄/q`

The data are the `κ(p/1009)` values for `p ≤ q/2`, extended to all units by `κ(q − p) = conj κ(p)`. Each value is placed at `x̃ = p̄/q`, and its Fourier coefficients are compared with those of the same values placed at `p/q`.

| `m` | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| `\|c_m\|` in `p̄/q` | `1.88e-2` | `1.30e-2` | `5.3e-3` | `5.9e-3` | `1.0e-3` | `4.7e-3` |
| `\|c_m\|` in `p/q` | `1.8e-4` | `1.5e-3` | `2.2e-4` | `1.8e-3` | `3.4e-4` | `3.1e-4` |

This is the V18 regression (`R² = 0.9984` on `x*`) seen from the other side. In the variable `p/q`,
`κ` is spectrally almost flat. The small `m = 2, 4` modes in `p/q` are the bounded-`p` bulbs
(`p = 2` sits at `x̃ ≈ ½`).

## B2 — Ramanujan sums in the spectrum, and the jump law at every denominator

**Observation.** Write `Im κ(x̃) ≈ Σ S_m sin(2πm x̃)`. The products `m·S_m` do not decay smoothly. They are
large at highly composite `m` and suppressed at primes ≥ 5:

| `m` | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 11 | 12 | 60 |
|---|---|---|---|---|---|---|---|---|---|---|
| `m·S_m` | `−.022` | `−.033` | `−.023` | `−.034` | `−.012` | `−.042` | `−.007` | `−.006` | `−.054` | `−.076` |

**Mechanism.** Let a function have a jump `J(p'/q')` at every rational and be smooth otherwise. Each jump contributes
`−J·((x − p'/q'))`, and `((y)) = −(1/π) Σ sin(2πmy)/m`. Summing over the reduced `p'`, with jumps
depending on `q'` only, gives

    π m S_m = Σ_{q'} J_{q'} c_{q'}(m) + (smooth part),        c_{q'}(m) = Ramanujan sum.

At a prime `m ≥ 5`, `c_{q'}(m) = μ(q')` for every `q' < m`, so the contributions alternate in sign and cancel.
At composite `m` the terms with `q' | m` add up. That is the table above.

**Fit.** Least squares of `π m S_m` on `c_1(m), …, c_Q(m)` and a `1/m²` background. The bounded-`p` bulbs (`p ≤ 4`) are excluded, since they are C1‴'s separate term.

| `m`-window | `Q` | R² | `J·q'²`, `q' = 1 … 8` |
|---|---|---|---|
| [4, 60] | 10 | 0.976 | `−.074 −.146 −.199 −.247 −.241 −.212 −.177 −.134` |
| [8, 120] | 12 | 0.926 | `−.064 −.129 −.196 −.250 −.273 −.255 −.299 −.260` |
| [12, 200] | 14 | 0.763 | `−.039 −.112 −.201 −.259 −.286 −.267 −.322 −.301` |
| [20, 300] | 16 | 0.546 | `+.009 −.073 −.205 −.289 −.285 −.284 −.347 −.357` |
| V21 (direct windows) | | | `—  −.136 −.167 −.168 −.122 −.104 (−.020)` |

**Reading.**

1. For `q' = 3 … 6`, `J·q'² ∈ [−0.29, −0.20]` in every window. The jump law is `q'^{-2}`-type, as V21
   said, with coefficient about `−0.25`. There may be a slow extra factor: `J·q'²` drifts upward in size with
   `q'` in the wider windows.
2. V21's direct values are smaller in size (`−0.10` to `−0.17`). This is expected, because V21 averages
   over windows `0.002 < |δ| < 0.009` that cut into the linear approach. V21's `q' = 7` value (`−0.0004`)
   is the "window-limited" artefact its own caveat names. Spectrally, `q' = 7` shows no collapse.
3. `q' = 1, 2` are least stable, because they are the most confounded with the smooth background.
4. R² falls at high `m`. At a single `q = 1009` the modes above `m ≈ 100` are noise-dominated.
   Only the low windows carry weight.

**Negative half.** The same model for `Re κ` uses the basis `B̄₂(x − p'/q')`, whose Fourier
coefficients are `∝ c_{q'}(m)/m²`, so a V-cusp has a slope jump `D_{q'}`. It fits poorly
(R² = 0.53, 0.34, 0.32, 0.36), and `m²C_m` grows (`−0.41` at `m = 60`). So the
rational singularities of `Re κ` are **not** V-cusps: they are softer than a jump and harder than a V.
~~The profile `|δ| log|δ|`, the Hilbert partner of a jump, is the natural next candidate.~~ *Update 2026-09-29:* the harmonic conjugate of a jump is `log|δ|`, and at `q = 2003` the profile is log type, with Brjuno-type organisation (register V33). C14″ currently says "V-cusps"; that wording
is not supported at this resolution.

**Why it matters.** Divisor and Ramanujan arithmetic in Fourier coefficients is the
signature of Eisenstein-type objects. Examples: Davenport's `Σ ((nx))/n^s` has sine coefficients
`σ_{1−s}(m)/m`, and weight-2 Eisenstein series have `σ₁`. B2 turns C14″'s "Fourier form with Ford weights" into
a computable identity, `πmS_m = Σ J_{q'} c_{q'}(m)`. Proving C14″ reduces to showing the jump law
`J_{q'}`, one number per denominator, and next move 2 of the register becomes a spectral fit
at `q = 2003` and beyond.

## B3 — the wake denominator `2^q − 1` divides the norm of the parabolic coefficient

At `λ = ζ_q`, `f^q(z) = z + a_q z^{q+1} + …` with `a_q = a_q(ζ_q)`, where `a_q(λ) ∈ ℤ[λ]`, and
`ι = b/a_q²` in normal form. So the primes of `N(a_q) = Res(Φ_q, a_q)` govern the
denominators of `ι`. V13's `ι_{1/3}`, `ι_{1/4}` and `ι_{1/5}` have denominators `21²`, `68²` and `5²·11²·31`.

**Lemma (proven for `λ ∈ {2, 4}`).** Let `ℓ > 2q + 2` be prime and `λ ∈ {2, 4}` with
`ord_ℓ(λ) = q`. Then `ℓ | N(a_q)`.

*Proof.*
1. `ℓ ≡ 1 (mod q)` splits completely in `ℚ(ζ_q)`. The ring maps `ℤ[ζ_q] → 𝔽_ℓ` send `ζ_q` to the
   primitive `q`-th roots of unity of `𝔽_ℓ`, and `λ` is one of them. Pick the prime `𝔩` with `ζ_q ≡ λ (mod 𝔩)`.
   Then `a_q(ζ_q) ≡ a_q(λ) (mod 𝔩)`.
2. **`λ = 2`.** `2z + z² = (1 + z)² − 1`, so `f^q(z) = (1 + z)^{2^q} − 1` and `a_q(2) = C(2^q, q + 1)`.
   Its numerator `2^q (2^q − 1) ⋯` contains `2^q − 1 ≡ 0 (mod ℓ)`, and `ℓ ∤ (q + 1)!`.
3. **`λ = 4`.** `4z + z²` is conjugate by `w = z + 2` to the Chebyshev map `w² − 2`, at its fixed point
   `w = 2`. With `N = 2^q`, `f^q` corresponds to `2T_N(w/2)`, and
   `2T_N(1 + z/2) = Σ_k 2N (N + k − 1)!/((N − k)!(2k)!) z^k`. At `k = q + 1`, the product
   `(N + k − 1)!/(N − k)!` contains `(N − 1)(N + 1) = 4^q − 1 ≡ 0 (mod ℓ)`, and `ℓ ∤ (2q + 2)!`.
4. So `a_q(ζ_q) ∈ 𝔩` and `ℓ | N(a_q)`. ∎

Both closed forms are checked as exact integer identities: `a_q(2)` for `q ≤ 24`, `a_q(4)` for `q ≤ 11`.

**Validated extension.** `λ = −2` is conjugate to the same Chebyshev map, but at its interior fixed
point `w = −1 = 2cos(2π/3)`, and the endpoint formula does not apply. Empirically, for every `q ≤ 20` and every prime
`ℓ > 2q + 2` dividing `Φ_q(2)Φ_q(4)Φ_q(−2)`:

    v_ℓ(N(a_q)) = #{λ ∈ {2, 4, −2} distinct mod ℓ : ord_ℓ(λ) = q}.

This gives **27 checks out of 27, all with equality.** For example, `31², 127², 8191²` come from `2` and `4` (Mersenne primes, `q` odd).
`43², 683², 2731²` come from `4` and `−2` (`ord_ℓ 2 = 2q`), and `257`, `65537` from `4` alone.

**Why it matters.** `M = 2^q − 1` is the denominator of the external angles: the `p/q` wake has
width `1/M` (P4), and P7–P9 work mod `M`. B3 shows that the same number sits in the cyclotomic
arithmetic of the analytic invariant `ι`. The mechanism is reduction mod `𝔩`: the parabolic multiplier
`ζ_q` becomes `2` or `4`, and the parabolic germ becomes a repelling fixed point of an exactly
solvable map, the power map or a Chebyshev map.

This is the first arithmetic statement about `ι` beyond "`ι ∈ ℚ(ζ_q)`". It also suggests a
**dynamics-over-𝔽_ℓ bridge** to `finite-mandelbrot-research`. There, `f_c` mod `ℓ` at
`c ≡ c_root(p/q)` could be read at the primes where `ζ_q ≡ 2`.

**Open.**
- A proof for `λ = −2`.
- Why equality holds, with no extra `ℓ`.
- The sporadic primes of `N(a_q)`, for example `987211` at `q = 7` and `12073` at `q = 8`, which none of the three multipliers explains.

## B4 — Dedekind sums are not the two-ended variable (negative)

Rademacher's `12 s(p, q) ≈ a₁ − a₂ + …` sees both ends of the continued fraction, so `s(p, q)` was
the obvious candidate for C16′'s "two resonance ends". It fails. The residual of `G` after a
100-bin function of `x*` has sd `1.63e-3`:

| regressor | residual sd | variance removed |
|---|---|---|
| none | `1.628e-3` | — |
| `1/p` | `1.618e-3` | 1.2% |
| `s(p,q)/q` | `1.595e-3` | 3.9% |
| both | `1.588e-3` | 4.9% |

`corr(residual, s) = +0.20`, so `r² = 3.9%`. With the `p ≤ 4` bulbs removed it is `+0.06`, which is 0.4% of the variance.
So `s(p, q)` does better than `1/p`, but most of its signal is the bounded-`p` bulbs, and over 95% of the
residual is untouched.

`s(p, q)` is the vendored `rational_dynamics_py.dedekind_sum`. The local `dedekind` this spike first used
skipped dividing out `gcd(h, k)` and was wrong whenever `gcd(h, k) > 1` (`s(2, 4) = −1/32`, not `0`); with the
prime `q = 1009` every `gcd(p, q) = 1`, so the numbers above are unaffected.

---

## 5. Ranking

1. **B3** is the most useful bridge. It is proven, it is new, and it ties the two halves of the study (rays mod `M`,
   and `ι` in `ℚ(ζ_q)`) with a one-paragraph proof. It opens a finite-field route.
2. **B2** is next. It turns the open cusp/jump laws into a spectral computation. It already gives the jump law at
   `q' = 3 … 8` in one fit, and it shows that C14″'s "V-cusp" wording needs revision.
3. **B1** confirms the variable. It is cheap, decisive, and belongs next to V18.
4. **B4** is closed as a negative. It should not be retried.

## 6. Register entries (applied 2026-09-28)

- **P13** (the B3 lemma, for `λ ∈ {2, 4}`), with `kernel/bulbford/norms.py` and `tests/test_parabolic_norms.py`.
- **V31** (B1 and B2 at `q = 1009`: the `p̄`-spectrum, the Ramanujan-sum jump law, and the failure of the V-cusp model for `Re κ`).
- **V32** (the `λ = −2` extension and the equality `v_ℓ = #zeros`, 27/27 for `q ≤ 20`).
- **C14″:** replace "V-cusps" by "rational singularities of `Re κ` (profile open; not V-type at `q = 1009`)".
- **Next move 2:** redo B2 at `q = 2003`, and test the `|δ| log|δ|` profile for `Re κ`. *Done 2026-09-29 (V33): log type, not `|δ| log|δ|`.*
