"""Integer-relation search for the horn-map constants (paper §9): each of Re x, Im x, |x|², arg x/π against every
subset of ≤ 3 of {1, π, π², 1/π, log 2, ζ(3), Catalan, γ}; PSLQ at 28 digits, coefficients ≤ 1000; a hit is
re-checked at 32 digits.  Writes data/pslq.txt."""
from itertools import combinations
import mpmath as mp

mp.mp.dps = 34
CONSTANTS = {
    "kappa_0": mp.mpc("0.023825887402200569000729189645402", "-0.052304659114003666316522486207571"),
    "kappa_0^(3)": mp.mpc("0.143599278095291569488064104088", "-0.030311647898502537611943110858"),
    "K_1/2": mp.mpc("0.018405261617964781008961965439681", "-0.042903723880230475207353337798906"),
}
BASIS = {"1": mp.mpf(1), "pi": mp.pi, "pi^2": mp.pi ** 2, "1/pi": 1 / mp.pi, "log2": mp.log(2),
         "zeta3": mp.zeta(3), "Catalan": mp.catalan, "gamma": mp.euler}
with open("data/pslq.txt", "w") as out:
    for name, x in CONSTANTS.items():
        found = []
        for qname, v in (("Re", mp.re(x)), ("Im", mp.im(x)), ("|.|^2", abs(x) ** 2), ("arg/pi", mp.arg(x) / mp.pi)):
            for k in range(1, 4):
                for sub in combinations(BASIS, k):
                    vec = [v] + [BASIS[s] for s in sub]
                    with mp.workdps(28):
                        rel = mp.pslq(vec, maxcoeff=1000, maxsteps=20000)
                    if rel and abs(mp.fsum(c * t for c, t in zip(rel, vec))) < mp.mpf(10) ** -30:
                        found.append((qname, sub, rel))
        line = f"{name}: {len(found)} relations" + (f" {found}" if found else "")
        out.write(line + "\n"); print(line, flush=True)
