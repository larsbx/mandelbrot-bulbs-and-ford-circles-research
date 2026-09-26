# Exact versus numerical evidence boundary

The research register mixes several evidence grades deliberately: proved
identities, rigorous ball computations, floating-point continuation/FFT
experiments, conjectures, and falsified hypotheses.  This note fixes the
boundary for the new cross-repository exact replay work.

## Exact low-q corpus

`data/parabolic_index_exact_vectors.json` pins four exact Gaussian-rational
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

- cycle continuation in `bulbford/dynamics.py`;
- Newton solves for centers and multiplier targets;
- Cauchy/FFT coefficient extraction in `bulbford/taylor.py`;
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
are copied byte-for-byte to `data/cyclotomic_germ_v1_vectors.json`
(`q <= 8`, every unit `p`).

`tests/test_cyclotomic_replay.py` replays them against two local instruments
that share no code with the kernel: `scripts/exact_index.py` (exact equality
of every coordinate) and the Arb ball of `bulbford/index.py` (containment).
Galois equivariance `ι_{p/q} = σ_p(ι_{1/q})` is an exact test in the kernel.

The kernel stage is a reference, not yet the canonical Mojo kernel, and it
does not name its coefficient.  Reading `[w^q] 1/P` as the holomorphic index
is this repository's.  For `q > 8`, `bulbford/index.py` remains a rigorous
ball computation rather than a canonical exact coefficient object.

## Joint box certificates for bulb centres

`bulbford/certify.py` applies the certificate calculus of
`larsbx/finite-mandelbrot-research` (`docs/finite-certificate-calculus.md`)
to satellite centres.  A centre of period `q` has critical-orbit type
`(0, q)`: `R_{0,q} = Q_q` is the Gleason polynomial, and a certificate is one
dyadic complex box `β` with

```text
K_{Q_q}(β) ⊂ int β                                      (unique simple root)
0 ∉ Q_j(β) − Q_i(β)   for all (i,j) ∈ F_{0,q}(q)         (exact period q)
```

checked on the same `β` in rational interval arithmetic with outward dyadic
rounding.  `data/center_certificates.json` stores 79 accepted boxes
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

## Wake combinatorics

`bulbford/wake.py` builds, in exact `Q/Z` arithmetic, the doubling cycle of
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
`data/rational_dynamics_r1_vectors.json`. The consumer replay covers
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
