"""Iterate the renormalization hierarchy G₁ = 𝒫₀ = 𝒫(z + z²), G_{k+1} = 𝒫(conj G_k), and report
κ_k = ι(G_k) − ½ = a₂/(2πi a₁²) of the horn map producing G_k.  Converges (if 𝒫 is contracting) to the index of the
parabolic-renormalization fixed point.  Writes data/levels.tsv (k, Re κ_k, Im κ_k, horn a₁, time)."""
import pickle, time, sys
import mpmath as mp
from bulbford.germ import horn_germ
from bulbford.horn import horn_coeffs
DPS, M, D = 60, 24, 40
mp.mp.dps = DPS
a = pickle.load(open("data/horn_coeffs_z2_M30_dps90.pkl", "rb"))[:M + 1]
out = open("data/levels.tsv", "w")
def record(k, a):
    with mp.workdps(DPS):
        kap = a[2] / (2j * mp.pi * a[1] ** 2)
    line = f"{k}\t{mp.nstr(kap.real, 30)}\t{mp.nstr(kap.imag, 30)}\t{mp.nstr(a[1], 15)}"
    out.write(line + "\n"); out.flush(); print(line, flush=True)
record(1, a)
for k in range(2, int(sys.argv[1]) + 1 if len(sys.argv) > 1 else 7):
    t0 = time.time()
    T = horn_germ(a, D, 1)                                  # Taylor coefficients of G_{k−1}
    A = T[2]
    germ = tuple(mp.conj(T[j] / A ** (j - 1)) for j in range(2, D + 1))   # normalized, conjugated
    a, al = horn_coeffs(M, h=0.5, N=64, R=600, n=800, dps=DPS, K=36, germ=germ)
    a = [mp.mpc(x) for x in a]
    record(k, a)
    print(f"   level {k}: {time.time()-t0:.0f}s, aliasing {mp.nstr(al, 2)}, |a_M| = {mp.nstr(abs(a[M]), 3)}", flush=True)
