"""Iterate the renormalization hierarchy G₁ = 𝒫₀ = 𝒫(z + z²), G_{k+1} = 𝒫(conj G_k); record κ_k = ι(G_k) − ½ =
a₂/(2πi a₁²).  Usage: levels.py TAG h M dps K_levels [holo].  With `holo` the conjugation is dropped:
G_{k+1} = 𝒫(G_k), the holomorphic operator of Inou–Shishikura / Lanford–Yampolsky (fixed point f_*).  Higher levels have germs with small discs of convergence
(level 2: |v| < 2.6 in normalized units), so the horn map must be sampled high (h ≈ 1: orbits stay at |v| ≲ 1);
the noise e^{2πnh}·10^{−dps} in a_n is absorbed by dps.  Writes data/levels_<TAG>.tsv and per-level coefficients."""
import pickle, time, sys
import mpmath as mp
from bulbford.germ import horn_germ
from bulbford.horn import horn_coeffs
TAG, h, M, DPS, KMAX = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
twist = (lambda x: x) if sys.argv[6:] == ["holo"] else mp.conj
D = 40
mp.mp.dps = DPS
a = pickle.load(open("data/horn_coeffs_z2_M30_dps90.pkl", "rb"))[:M + 1]
out = open(f"data/levels_{TAG}.tsv", "w")
def record(k, a, extra=""):
    with mp.workdps(DPS):
        kap = a[2] / (2j * mp.pi * a[1] ** 2)
        T = horn_germ(a, D, 1); A = T[2]
        g = [abs(T[j] / A ** (j - 1)) for j in range(2, D + 1)]
        rad = 1 / max(x ** (1 / (j + 1)) for j, x in enumerate(g) if j >= 15 and x > 0)   # radius estimate in v
    line = f"{k}\t{mp.nstr(kap.real, 30)}\t{mp.nstr(kap.imag, 30)}\t{mp.nstr(a[1], 15)}\t{mp.nstr(rad, 3)}\t{extra}"
    out.write(line + "\n"); out.flush(); print(line, flush=True)
    pickle.dump(a, open(f"data/levels_{TAG}_a{k}.pkl", "wb"))
record(1, a)
for k in range(2, KMAX + 1):
    t0 = time.time()
    T = horn_germ(a, D, 1); A = T[2]
    germ = tuple(twist(T[j] / A ** (j - 1)) for j in range(2, D + 1))   # normalized, conjugated unless holo
    a, al = horn_coeffs(M, h=h, N=64, R=600, n=800, dps=DPS, K=36, germ=germ)
    a = [mp.mpc(x) for x in a]
    record(k, a, f"{time.time()-t0:.0f}s alias {mp.nstr(al, 2)} |a_M| {mp.nstr(abs(a[M]), 3)}")
