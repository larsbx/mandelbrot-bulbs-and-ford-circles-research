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
coefficient representation.  Until that shared representation is implemented
and independently replayed, general-q values from `bulbford/index.py` remain
rigorous ball computations rather than canonical exact coefficient objects.
