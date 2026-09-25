"""Level-2 test, route 2: horn map of the germ 𝒫₀ (and of its complex conjugate), i.e. parabolic renormalization
applied twice to z + z².  𝒫₀(W) = W + A W² + …, A = 2πi a₁; normalized germ g(v) = A𝒫₀(v/A) = v + v² + Σ (T_j/A^{j−1}) v^j."""
import pickle, time
import mpmath as mp
from bulbford.germ import horn_germ
from bulbford.horn import horn_coeffs
mp.mp.dps = 60
a = pickle.load(open("data/horn_coeffs_z2_M30_dps90.pkl", "rb"))
D = 40
T = horn_germ(a, D, 1)                                    # Taylor coefficients of 𝒫₀ up to degree D
A = T[2]
germ = tuple(T[j] / A ** (j - 1) for j in range(2, D + 1))      # (1, g₃, g₄, …)
germ_c = tuple(mp.conj(x) for x in germ)                  # 𝒫₀ conjugated by z ↦ z̄
print("g₃ =", mp.nstr(germ[1], 12), " résit(𝒫₀) = 1 − g₃ ⇒ ι(𝒫₀) = g₃ =", mp.nstr(germ[1], 12), "(κ₀ + ½ = 0.5238… − 0.0523…i)")
for name, gm in (("P0", germ), ("conj P0", germ_c)):
    for h in (0.5, 0.75):
        t0 = time.time()
        b, al = horn_coeffs(3, h=h, N=64, R=600, n=800, dps=50, K=36, germ=gm)
        with mp.workdps(50):
            k2 = b[2] / (2j * mp.pi * b[1] ** 2)
            print(f"{name:>8} h={h}: a0={mp.nstr(b[0], 3)}  a1={mp.nstr(b[1], 14)}  a2={mp.nstr(b[2], 14)}  "
                  f"a2/(2πi a1²) = {mp.nstr(k2, 16)}   conj = {mp.nstr(mp.conj(k2), 16)}  aliasing {mp.nstr(al, 2)}  {time.time()-t0:.0f}s", flush=True)
