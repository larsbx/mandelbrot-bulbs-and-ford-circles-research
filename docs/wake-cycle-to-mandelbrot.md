# Wake cycle to the Mandelbrot set

This note makes the parameter-plane interpretation of `kernel/bulbford/wake.py`
visual without changing the repository's evidence boundary.

`wake.py` is a thin adapter over the vendored `rational_dynamics_py`
(`vendor/python`, from `larsbx/finite-math-kernels`, pinned in `vendored.toml`);
it keeps the names `wake`, `mechanical`, `rotation_cycle` and `double` that
vizops reads, and puts this checkout's `vendor/python` on `sys.path` when it is
executed on its own.

The visuals are drawn in
[`larsbx/math-vizops`](https://github.com/larsbx/math-vizops), which reads this
repository and writes nothing back:

- a static 3/7 still,
  [`wiki/images/wake-cycle-3-7.svg`](https://github.com/larsbx/math-vizops/blob/main/wiki/images/wake-cycle-3-7.svg),
  whose printed angles vizops's estate test holds to `kernel/bulbford/wake.py`;
- an interactive page, `python -m vizops page wake-to-mandelbrot`, which
  renders a numerical Mandelbrot raster for `z² + c`, lets you choose any
  reduced `p/q` with `2 ≤ q ≤ 12`, and shows the exact finite angle data
  computed by `kernel/bulbford/wake.py` itself, embedded at build time.

## What maps to what

The map is **not**

```text
q points of the doubling cycle  ->  q parameter points on the bulb
```

Instead the repository's construction is:

```text
rotation number p/q
        |
        v
exact q-cycle of θ -> 2θ mod 1
        |
        | choose the shortest adjacent arc
        v
characteristic pair (θ₋, θ₊)
        |
        | imported rational parameter-ray landing [DH/Mil00]
        v
two external parameter rays
        |
        v
one common landing point = root of B_{p/q}
        |
        v
region between the rays = the p/q wake
```

The other `q - 2` angles are still essential: they determine the finite
rotation orbit and therefore which adjacent pair is characteristic. They are
not additional bulb-root landing points.

## Exact specimen: 3/7

For `p/q = 3/7`, `kernel/bulbford/wake.py` gives the doubling orbit with
denominator `2⁷ - 1 = 127`.

In orbit order:

```text
21 -> 42 -> 84 -> 41 -> 82 -> 37 -> 74 -> 21    (mod 127)
```

In angular order:

```text
21, 37, 41, 42, 74, 82, 84
```

Hence the characteristic pair is exactly

```text
θ₋ = 41/127
θ₊ = 42/127
width = 1/127.
```

For the standard quadratic parameter plane, the main-cardioid boundary point
with multiplier

```text
λ = exp(2πi · 3/7)
```

has

```text
c = λ/2 - λ²/4
  ≈ -0.6063568844 + 0.4123997402 i.
```

The visual marks this numerical coordinate as the attachment root of the
`3/7` satellite. The coordinate is only a display aid here; it is not used
to promote any numerical observation to an exact claim.

## Evidence boundary

The visuals deliberately use three visual grades:

| grade | shown in the visual | status |
|---|---|---|
| exact finite arithmetic | rotation word, doubling orbit, `θ₋`, `θ₊`, characteristic width | computed by `kernel/bulbford/wake.py` and embedded by vizops |
| imported dynamics | the two rational parameter rays land together at the root of `B_{p/q}` | `[DH/Mil00]`, already named by the kernel |
| visual aid only | Mandelbrot raster, displayed root coordinate, dashed ray traces, satellite outline | numerical or schematic; no certificate claim |

In particular, the dashed curves are **not** numerical integrations of
parameter rays. Their only purpose is to make the combinatorial-to-parameter
mapping visually obvious.

## Why this matters for the Ford-circle comparison

Ford circles already make the `p/q` hierarchy visible on the rational line:
Farey adjacency is tangency and the radius is exactly `1/(2q²)`.

The wake construction supplies the corresponding exact combinatorial address
on the Mandelbrot side. The finite research program can therefore keep the
layers separate:

```text
Farey / Ford arithmetic
        |
        +---- exact p/q combinatorics
        |
doubling wake arithmetic
        |
        +---- exact characteristic angles
        |
imported landing theorem
        |
Mandelbrot bulb root and wake
        |
finite bulb geometry / index corrections
```

This separation is useful because the Ford-scale analogy can be visualized
without treating the stronger geometric correction claims as already proved.
