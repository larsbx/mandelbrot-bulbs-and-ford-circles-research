"""V47: the first-order correction in 1/q is kinematic.  R_{p/q}(u) (bulb Taylor data, Cauchy on |u| = 1.5) against
ρ_{p,r}(ũ), ũ = u/(1 + u/(2πipq)) (renormalized family); residuals of r₂…r₈ and their ratio per doubling of N.
Degree 2 on |u| = 1.5, degree 3 (cubic cusp germ) on |u| = 1 (E-h).  Writes data/kinematic.txt."""
import numpy as np
from bulbford.taylor import taylor
from bulbford.renorm import unfolding_taylor, lavaurs_compose
from bulbford.dynamics import MAIN2, MAIN3
from bulbford.horn import QUADRATIC, CUBIC

K = 9
with open("data/kinematic.txt", "w") as out:
    cases = [(MAIN2, QUADRATIC, 1.5, p, r) for p, r in ((1, 1), (2, 1), (3, 1), (3, 2))] + [(MAIN3, CUBIC, 1.0, 1, 1)]
    for fam, germ, rad, p, r in cases:
        rho = unfolding_taylor(p, r, radius=rad, germ=germ)[:K]
        out.write(f"degree {fam.d}, p = {p}, r = {r}: |r_k − [u^k]ρ(ũ)|, k = 2…{K - 1}  (plain |r_k − ρ_k| for k = 2, 3 in brackets)\n")
        prev = None
        for N in (64, 128, 256, 512, 1024):
            q = p * N + r
            c = taylor(p, q, fam, r=rad, N=128).coeffs[:K]
            res = np.abs(c - np.array(lavaurs_compose(rho, p, q)))[2:]
            naive = np.abs(c - rho)[2:4]
            line = f"  q = {q:5d}: " + " ".join(f"{x:.2e}" for x in res) + "  [" + " ".join(f"{x:.2e}" for x in naive) + "]"
            if prev is not None:
                line += "   ratio " + " ".join(f"{a / b:.2f}" for a, b in zip(prev, res))
            out.write(line + "\n"); print(line, flush=True)
            prev = res
