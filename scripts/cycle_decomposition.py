"""ι_{p/q} = −1/(1−(2−λ₀)^q) − q Σ_C 1/(1−ρ_C) over the other q-cycles (index sum [Mil]).
All fixed points of f^q (f = λ₀z+z²) from Newton with many starts, deduplicated; sorted contributions."""
import sys, cmath
import numpy as np
from bulbford.index import index_complex

def fixed_points(lam, q, seeds=80000, seed=0):
    rng = np.random.default_rng(seed)
    z = (rng.uniform(-2, 2, seeds) + 1j * rng.uniform(-2, 2, seeds)).astype(complex)
    for _ in range(300):
        w, A = z.copy(), np.ones_like(z)
        for _ in range(q):
            A = (lam + 2 * w) * A; w = lam * w + w * w
        step = (w - z) / (A - 1); z = z - np.where(np.abs(step) < 1, step, step / np.abs(step))
    w, A = z.copy(), np.ones_like(z)
    for _ in range(q):
        A = (lam + 2 * w) * A; w = lam * w + w * w
    ok = np.abs(w - z) < 1e-10
    pts = []
    for zz, aa in zip(z[ok], A[ok]):
        if all(abs(zz - u) > 1e-7 for u, _ in pts): pts.append((zz, aa))
    return pts

for pq in sys.argv[1:] or ["1/7", "2/7", "1/11", "3/11"]:
    p, q = map(int, pq.split("/")); lam = cmath.exp(2j * cmath.pi * p / q)
    pts = fixed_points(lam, q)
    far = [(z, A) for z, A in pts if abs(A - 1) > 1e-3 and abs(z - (1 - lam)) > 1e-9]   # exclude 0 and β
    beta = 1 / (1 - (2 - lam) ** q)
    # group into cycles by multiplier value
    cyc = {}
    for z, A in far: cyc.setdefault(round(A.real, 6) + 1j * round(A.imag, 6), []).append(z)
    contrib = sorted(((-q / (1 - A), A, len(zs)) for A, zs in cyc.items()), key=lambda c: -abs(c[0]))
    total = -beta + sum(c[0] for c in contrib)
    print(f"\n{pq}: {len(pts)} fixed points found of 2^{q}={2**q} (incl. 0 counted once, β); {len(cyc)} other cycles "
          f"(expected {(2**q - 2 - q) // q}); ι(numeric sum)={total:.6f}  ι(series)={index_complex(p, q):.6f}  β-term={-beta:.2e}")
    for c, A, n in contrib[:8]:
        print(f"   |ρ|={abs(A):9.3f}  arg ρ/2π={cmath.phase(A) / (2 * cmath.pi):+.3f}  pts={n:>2}  −q/(1−ρ) = {c:+.5f}")
    print(f"   sum of remaining {max(0, len(contrib) - 8)}: {sum(c[0] for c in contrib[8:]):+.5f}")
