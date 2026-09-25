"""Asymptotic fit κ(1/q) = Σ_{j=0}^{J} c_j q^{-j} from data/p1_exact.tsv (mpmath, 60 digits).
Reports c_0 = κ₀ across (J, q_min) for stability, and d_j = c_j/(2πi)^j (real if κ(1/q) − κ₀ is a real-coefficient
series in ε = 2πi/q)."""
import sys
import mpmath as mp
mp.mp.dps = 60
SRC = sys.argv[1] if len(sys.argv) > 1 else "data/p1_nf.tsv"
rows = sorted((int(a), mp.mpc(mp.mpf(b), mp.mpf(c))) for a, b, c, *_ in (l.split("\t") for l in open(SRC) if l.strip()))
print(f"{len(rows)} values, q = {[q for q, _ in rows]}")
def fit(qmin, J):
    R = [(q, k) for q, k in rows if q >= qmin]
    if len(R) < J + 2: return None
    A = mp.matrix([[mp.mpf(q) ** (-j) for j in range(J + 1)] for q, _ in R]); y = mp.matrix([k for _, k in R])
    c = mp.lu_solve(A.T * A, A.T * y)
    res = max(abs((A * c)[i] - y[i]) for i in range(len(R)))
    return c, res, len(R)
print("\n q_min  J  n |              Re κ₀                    Im κ₀              | max resid")
best = {}
for qmin in (64, 128, 192, 256, 384):
    for J in (4, 6, 8, 10):
        f = fit(qmin, J)
        if not f: continue
        c, res, n = f; best[(qmin, J)] = c
        print(f" {qmin:>5} {J:>2} {n:>2} | {mp.nstr(c[0].real, 22):>24} {mp.nstr(c[0].imag, 22):>24} | {mp.nstr(res, 3)}")
c = best.get((128, 8)) or next(iter(best.values()))
print("\nq_min=128, J=8:  j | c_j | d_j = c_j/(2πi)^j")
for j in range(len(c)):
    d = c[j] / (2j * mp.pi) ** j
    print(f"  {j} | {mp.nstr(c[j], 12)} | {mp.nstr(d, 12)}")
