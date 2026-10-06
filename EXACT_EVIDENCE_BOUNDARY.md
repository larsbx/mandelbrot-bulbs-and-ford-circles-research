# Exact versus numerical evidence boundary

The research register mixes several evidence grades deliberately: proved
identities, rigorous ball computations, floating-point continuation/FFT
experiments, conjectures, and falsified hypotheses.  This note fixes the
boundary for the new cross-repository exact replay work.

## Exact low-q corpus

`tests/vectors/parabolic_index_exact_vectors.json` pins four exact Gaussian-rational
values of the parabolic fixed-point index:

```text
q=1, p=1  -> 0
q=2, p=1  -> 1/8
q=4, p=1  -> (1447 - 365 i)/4624
q=4, p=3  -> (1447 + 365 i)/4624
```

The file contains rational numerators and denominators only.  It contains no
floating-point coordinates and no tolerance.

These fixtures are intended for independent exact replay by finite-algebra
implementations.  The repository's existing Arb calculation may agree with
them, but that agreement is a cross-check rather than the reason the exact
fractions are exact.

## Numerical research lane

The following remain research instruments rather than exact certificate
objects:

- cycle continuation in `kernel/bulbford/dynamics.py`;
- Newton solves for centers and multiplier targets;
- Cauchy/FFT coefficient extraction in `kernel/bulbford/taylor.py`;
- dense asymptotic sweeps and fitted limits.

Their output remains VALIDATED numerical evidence under the research register.
No exact replay promotes those observations automatically.

## Exact algebra lane

The exact lane may establish finite statements such as:

- a root-of-unity relation in an exact coefficient ring;
- a finite truncated iterate;
- a vanishing order;
- a reciprocal-series coefficient when the constant coefficient is
  invertible;
- Galois/conjugation identities on exact coefficients.

It does not by itself establish bulb diameter, global component geometry,
asymptotic convergence, or any classical landing/connectivity theorem.

## General-q frontier

The q=1,2,4 vectors live in Q(i).  General p/q requires an exact cyclotomic
coefficient representation.  `larsbx/finite-math-kernels` now carries an
independent reference for it (bridge stages C1, C2, Q1, Q2 in
`tools/cyclotomic_reference.py`, commit `16bf936`), and its pinned vectors
are copied byte-for-byte to `tests/vectors/cyclotomic_germ_v1_vectors.json`
(`q <= 8`, every unit `p`).

`tests/test_cyclotomic_replay.py` replays them against two local instruments
that share no code with the kernel: `experiments/scripts/exact_index.py` (exact equality
of every coordinate) and the Arb ball of `kernel/bulbford/index.py` (containment).
Galois equivariance `ι_{p/q} = σ_p(ι_{1/q})` is an exact test in the kernel.

The kernel stage is a reference, not yet the canonical Mojo kernel, and it
does not name its coefficient.  Reading `[w^q] 1/P` as the holomorphic index
is this repository's.  For `q > 8`, `kernel/bulbford/index.py` remains a rigorous
ball computation rather than a canonical exact coefficient object.

## Joint box certificates for bulb centres

`kernel/bulbford/certify.py` applies the certificate calculus of
`larsbx/finite-mandelbrot-research` (`docs/finite-certificate-calculus.md`)
to satellite centres.  A centre of period `q` has critical-orbit type
`(0, q)`: `R_{0,q} = Q_q` is the Gleason polynomial, and a certificate is one
dyadic complex box `β` with

```text
K_{Q_q}(β) ⊂ int β                                      (unique simple root)
0 ∉ Q_j(β) − Q_i(β)   for all (i,j) ∈ F_{0,q}(q)         (exact period q)
```

checked on the same `β` in rational interval arithmetic with outward dyadic
rounding.  `experiments/data/center_certificates.json` stores 79 accepted boxes
(`2 <= q <= 16`, every unit `p`, half-width `2^-64`); the tests replay each
from its endpoints alone and check the boxes of one `q` are pairwise
disjoint.  The interval semantics is additionally replayed on the two boxes
pinned by the finite-mandelbrot reference plane (`c = -2`, 5/5 exclusions;
`M_{4,1}`, 18/18).

An unmet check is INCONCLUSIVE, never a disproof (the fail-closed direction
"refuse to conclude").  Seeds come from floating-point continuation and are
untrusted.

What the box does **not** carry: that the certified centre is the centre of
`B_{p/q}`.  That label rides on `[DH]` plus the continuation seed
(`SatelliteLabel`, VALIDATED), recorded as named tags in every record.

## Certified antipodes and G_ant enclosures

`kernel/bulbford/antipode.py` extends the joint box to the `ρ = −1` point. On a box
`Z × C` a two-variable Krawczyk inclusion for

```text
F(z, c) = (f_c^q(z) − z, (f_c^q)'(z) + 1)
```

gives one solution, and the type-`(0,q)` exclusions on the orbit of `Z`
give it exact period `q`. `ζ_q` is certified without an angle: all roots of
`X^q − 1` get disjoint Krawczyk boxes, and `ζ_q` is the upper-half-plane box
whose real part is strictly largest. With `λ₀ = ζ_q^p` and
`c_root = λ₀/2 − λ₀²/4` as interval expressions, the enclosure is taken for
the quadrance

```text
G_ant² = q⁴ Qd(c_ant − c_root) / Qd(1 − λ₀)
```

and `G_ant` is bracketed by rational square-root bounds.
`experiments/data/antipode_certificates.json` holds 79 accepted records (`q <= 16`),
each replayed from its endpoints; the stored bounds are recomputed, not
trusted.

`experiments/data/antipode_certificates_v14.json` extends this to every `p <= q/2` at
the V14 denominators `q = 59, 127, 251`. Along an orbit of length `q`,
rectangular arithmetic widens an enclosure by `|Re 2z_i| + |Im 2z_i|` per
step, not by `|2z_i|`, so the generator sizes each box from that growth
(estimated from the untrusted float orbit; it only chooses the box) and
stops any enclosure that passes `2^40` as INCONCLUSIVE. Only the `G_ant`
column of V14 becomes certified this way. The `|u_a|/2` column is certified
separately by `experiments/scripts/v14_certified.py`: on each certified `c_ant` box,
Arb encloses `λ_a = 1 − √(1 − 4c_ant)` (the branch `Re(1 − λ) > 0` is
checked on the ball, else refused) and `u_a = q² log(λ_a/λ₀)`. Arb's
square root and logarithm are rigorous balls, and `λ₀` is the angle-free
ball of the next section, so both V14
columns and `Δ` are certified; `experiments/data/v14_certified.json` holds the
brackets.

This is a certified enclosure of a finite algebraic quantity. Reading it as
the size proxy of `B_{p/q}` still needs `[DH]` and `SatelliteLabel`; it
says nothing about `q → ∞`, and so leaves every limit and conjecture in the
register at its current grade.

## Centre and antipode in one component

`kernel/bulbford/continuation.py` certifies, for all 79 bulbs with `q <= 16`, that
the P10 antipode lies on the boundary of the hyperbolic component of the P3
centre. Along the segment between the two certified boxes, a chain of
parameter boxes carries one attracting `q`-cycle: a parametric Krawczyk
inclusion on each box (one fixed point of `f_c^q` per `c`, analytic in
`c`), `|rho| < 1` on the interior boxes, junction inclusions so that
neighbouring boxes hold the same point, and on the last box a certified
positive derivative of `|rho|^2` along the path, so `|rho|` rises strictly to
`1` at the antipode. The computation uses Arb complex balls; the rational
certificate boxes enter exactly. `experiments/data/continuation_certificates.json`
holds the verdicts, `experiments/scripts/certify_continuation.py --check` recomputes
all of them, and CI recomputes `q <= 7`.

This removes the antipode half of `SatelliteLabel`. The centre half, that
the certified centre is the centre of `B_{p/q}`, remains a named import
outside the P14 range below.

## The centre's component touches the root

`kernel/bulbford/root.py` certifies, for every bulb with `q <= 4`, that the
hyperbolic component of the P3 centre has `c_root(p/q)` on its boundary.
`lambda0` comes from the polynomial/order interface `index.lambda_ball`.
Certified contour means evaluate argument-principle integrals around
`z0 = lambda0/2`. They count `q + 1` fixed points of `f_c^q` (and only
`alpha` for each proper divisor) in two nested disks, which isolates one
`q`-cycle `C(c)`. They also give its multiplier as a contour integral
that stays defined at the root, where `C` collides with `alpha`. Certified
signs of `d|rho|^2` and `d|lambda|^2` along the segment from the root to a
point `c1` show the segment lies in one hyperbolic component with the root
on its boundary. A continuation from the centre then reaches `C(c1)`.

The contour calculation also uses only algebraic roots. Write
`H(u) = u h(m + u)` for an integrand `h`. Boxes obtained by rotating a
rational rectangle with certified roots of unity cover the annulus
`r/a <= |u| <= r*a`, with rational `a > 1`. Finite ball evaluations there
exclude every denominator zero and principal-log cut and bound `|H|` by
`M`. Cauchy's Laurent coefficient bounds imply that the mean of `H` at
`r*zeta_N^k` differs from its constant coefficient, the normalized contour
integral, by at most `2M/(a^N - 1)`. That error is added to both coordinates
of the returned ball. Polynomial Taylor shifts reduce interval widening
near the multiple root. No angle, trigonometric function, or value of pi
is evaluated, including in the contour normalization. An uncertified
annulus or error bound is refused.

The log in the multiplier integral must be holomorphic on the whole inner
disk, not only on this annulus. The code verifies
`|2m/lambda0 - 1| + 2r' < 1` before any contour calculation and refuses otherwise.

What stays imported is one [DH] fact: only the main cardioid and `B_{p/q}`
have `c_root(p/q)` in their closure. The recertified coverage remains the
five fractions `1/2, 1/3, 2/3, 1/4, 3/4`. Each stored disk and endpoint is
replayed with its centre-to-endpoint continuation in the canonical tests;
`experiments/scripts/certify_root.py --check` also regenerates the corpus.
No `q >= 5` conclusion follows from this recertification. The previous
angle-based integration reported success at `2/5, 3/5` and refusals at
`1/5, 4/5` and `q >= 6`; those are historical observations about that
integrator, not results or limitations of the replacement method.

## The angle-free lane

No exact or certified computation in `kernel/` evaluates an angle, `π`, or a
trigonometric function. Each root of unity is built from its polynomial and
selected by order alone. There are three independent constructions, and the
tests make each pair that overlaps meet.

| Module | Field | Construction | What it gives |
|---|---|---|---|
| `kernel/bulbford/spread.py` | `Q` | `cos(πp/q)` is the `p`-th largest root of the integer Chebyshev polynomial `U_{q−1}` (since `sin qθ = sin θ · U_{q−1}(cos θ)`); a Sturm sequence over `Q` isolates it and bisection narrows it | rational brackets for `Qd(1 − λ₀) = 4(1 − c²)` and `pred² = Qd(1 − λ₀)/q⁴` |
| `kernel/bulbford/antipode.py` | `Q(i)` | Krawczyk boxes on `X^q − 1`; the untrusted float seeds now come from the `spread.py` brackets | `ζ_q` boxes |
| `kernel/bulbford/index.py` | Arb | Arb's certified roots of `Φ_q`, then the same order check | the `ζ_q` and `λ₀` balls behind the index and the `|u_a|/2` enclosures |

Notes on the table:

- **The rational route.** The spread polynomials of rational trigonometry
  come from the same polynomials: `S_n(1 − c²) = (1 − c²) U_{n−1}(c)²` and
  `S_n ∘ S_m = S_{nm}`. So the bulb–Ford first-order prediction is a
  bracketed algebraic number with no complex numbers at all. It is
  `Qd(1 − λ₀)` times the squared Ford diameter `(1/q²)²`.
- **Cross-checks.** `tests/test_spread.py` checks the `spread.py` brackets
  against the `Q(i)` boxes for `q <= 12`. `tests/test_angle_free_roots.py`
  checks the Arb balls against the `Q(i)` boxes for `q <= 16`.
- **P14.** `root.py` uses this same Arb root interface for `lambda0` and
  every contour node. It has no numerical-lane exemption. Its logarithm
  and exponential enclosures remain rigorous Arb computations; the
  no-angles rule does not turn these into exact rational values.

`tools/audit_angles.py` enforces the rule on Python tokens, so prose in
docstrings and comments stays free. The audit runs in CI.

- **Exempt modules.** The modules of the numerical research lane above are
  declared exempt, each with its reason: `dynamics.py`, `cycles.py`,
  `taylor.py`, `norms.py`, and `implosion.py`, where `π` belongs to the
  Fatou-coordinate statements themselves.
- **Stale exemptions fail.** A declared module that no longer uses an angle
  is a breach, so the list only shrinks deliberately.
- **P10/V14 certificates.** Regenerating these certificate files with the angle-free
  seeds changed no verdict, field or certificate box. 76 rational brackets
  changed their digits, all of them built through the `λ₀` box: 6 `G_ant`
  brackets in `antipode_certificates.json`, 18 in
  `antipode_certificates_v14.json`, and 52 in `v14_certified.json`. Each new
  bracket overlaps the old one.

## Wake combinatorics

`kernel/bulbford/wake.py` builds, in exact `Q/Z` arithmetic, the doubling cycle of
rotation number `p/q` and its characteristic arc `(θ₋, θ₊)`, of width
`1/(2^q − 1)`.  Finite facts, tested for `q <= 16`: the cycle acts as the
rotation `x_i ↦ x_{i+p}`; the conjugate wake is the mirror; wakes are
disjoint and ordered like Farey fractions; the mediant wake sits in the gap
between Farey neighbours.  Doubling agrees with the shared R1 vectors.
Uniqueness of the rotation cycle `[Gol92]` and landing of the two rays at the
root of `B_{p/q}` `[DH/Mil00]` are imports.


## Shared rational-dynamics contract

The arithmetic of the internal fraction is now cross-checked against
`larsbx/finite-math-kernels` R1 through
`tests/vectors/rational_dynamics_r1_vectors.json`. The consumer replay covers
reduction-compatible fractions, doubling modulo one, modular inverse, the
centered inverse representative, canonical continued fractions, convergent
denominators, and the specimen identity

```text
3/7 = [0;2,3],
signed inverse numerator = -2,
previous convergent denominator = 2.
```

This does not move the interpretation of `x*` or any bulb claim into the
shared kernel. The kernel owns only the finite rational arithmetic; this
repository owns the use of that arithmetic in its research register.

Since 2026-10 the finite arithmetic itself is the shared code, not a local
copy checked against it: `kernel/bulbford/{cf,wake,cycles}.py` are adapters
over the vendored `rational_dynamics_py` (`vendor/python`, pinned in
`vendored.toml`). The R1 vector replay above still runs, now through those
adapters; the exact/imported boundary of the wake combinatorics is unchanged.

Likewise the Krawczyk operator (Krawczyk 1969), its rounded-midpoint
preconditioner, the strict-interior Krawczyk-Moore test (Moore 1977) and box
disjointness are the vendored `root_isolation_py` over `closed_interval`
(finite-math-kernels `docs/root-isolation-spec.md`), run by
`kernel/bulbford/{certify,antipode}.py` on their own boxes, maps and
`2^-prec` outward rounding. Every stored certificate and every replayed
verdict is byte-identical to the local operator it replaced; what a passing
inclusion means, and the imports it needs, are unchanged.
