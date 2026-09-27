"""Order of limits at x̃ → 1/3⁻: G for x̃ = [0;3,N,a1] (p/q = [0;a1,N,3]) on a grid of (N, a1)."""
from paths import DATA
import json
from bulbford.cf import from_cf
from bulbford.taylor import kappa_fft
rows = []
for N in (32, 64, 128):
    for a1 in (2, 3, 5, 9, 16, 32, 64):
        p, q = from_cf((0, a1, N, 3))
        t = kappa_fft(p, q); G = abs(t.solve(-1, 2.0)) / 2
        rows.append(dict(N=N, a1=a1, p=p, q=q, kappa=[t.coeffs[2].real, t.coeffs[2].imag], G=G))
        print(f"N={N:>3} a1={a1:>2} p/q={p}/{q:<6} κ={t.coeffs[2]:+.5f} G={G:.5f}", flush=True)
    p, q = from_cf((0, N, 3)); t = kappa_fft(p, q); G = abs(t.solve(-1, 2.0)) / 2
    rows.append(dict(N=N, a1=None, p=p, q=q, kappa=[t.coeffs[2].real, t.coeffs[2].imag], G=G))
    print(f"N={N:>3} a1= ∞(p=3) p/q={p}/{q:<6} κ={t.coeffs[2]:+.5f} G={G:.5f}", flush=True)
json.dump(rows, open(DATA / "limit_grid.json", "w"))
