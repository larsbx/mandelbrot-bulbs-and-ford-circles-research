"""κ(1/q) = Σ_j c_j q^{-j}: free fits over (q_min, J), then the same with c₁ fixed to −i/(2π).
Reports κ₀ spread (stability) and the coefficient table of the best free fit."""
import sys, mpmath as mp
mp.mp.dps = 60
SRC = sys.argv[1] if len(sys.argv) > 1 else "data/p1_nf.tsv"
rows = sorted((int(a), mp.mpc(mp.mpf(b), mp.mpf(c))) for a, b, c, *_ in (l.split("\t") for l in open(SRC) if l.strip()))
C1 = -1j / (2 * mp.pi)

def fit(qmin, J, fix_c1=False):
    R = [(q, k - (C1 / q if fix_c1 else 0)) for q, k in rows if q >= qmin]
    js = [j for j in range(J + 1) if not (fix_c1 and j == 1)]
    if len(R) < len(js) + 2: return None
    A = mp.matrix([[mp.mpf(q) ** (-j) for j in js] for q, _ in R]); y = mp.matrix([k for _, k in R])
    c = mp.lu_solve(A.T * A, A.T * y); res = max(abs((A * c)[i] - y[i]) for i in range(len(R)))
    return dict(zip(js, c)), res, len(R)

print(f"{len(rows)} values, q = {rows[0][0]}…{rows[-1][0]}")
for fix in (False, True):
    print(f"\n{'c₁ fixed = −i/(2π)' if fix else 'free fit'}:   q_min  J  n | Re κ₀ | Im κ₀ | max resid")
    k0s = []
    for qmin in (64, 128, 256, 512):
        for J in (4, 6, 8, 10):
            f = fit(qmin, J, fix)
            if not f: continue
            c, res, n = f; k0s.append((res, c[0]))
            print(f"   {qmin:>5} {J:>2} {n:>2} | {mp.nstr(c[0].real, 20)} | {mp.nstr(c[0].imag, 20)} | {mp.nstr(res, 3)}")
    good = [k for r, k in k0s if r < 1e-20]
    if good:
        print(f"   κ₀ over fits with resid < 1e-20: Re ∈ [{mp.nstr(min(k.real for k in good), 18)}, {mp.nstr(max(k.real for k in good), 18)}], "
              f"Im ∈ [{mp.nstr(min(k.imag for k in good), 18)}, {mp.nstr(max(k.imag for k in good), 18)}]")
best = min((f for f in (fit(m, J) for m in (64, 128) for J in (8, 10)) if f), key=lambda f: f[1])
print("\nbest free fit coefficients:  j | c_j | c_j·(2π)^j/i^j")
for j, cj in sorted(best[0].items()):
    print(f"  {j} | {mp.nstr(cj, 15)} | {mp.nstr(cj * (2 * mp.pi) ** j / (1j) ** j, 15)}")
print(f"\nc₁ + i/(2π) = {mp.nstr(best[0][1] - C1, 5)}")
