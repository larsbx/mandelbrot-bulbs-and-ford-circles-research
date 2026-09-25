"""κ(1/q) for large q (C20 constant κ₀), FFT route. Writes data/p1_limit.json."""
import json, sys
from bulbford.taylor import kappa_fft
rows = []
for q in map(int, sys.argv[1:] or ["1009", "2003", "4001", "8009"]):
    t = kappa_fft(1, q); G = abs(t.solve(-1, 2.0)) / 2
    rows.append(dict(q=q, kappa=[t.coeffs[2].real, t.coeffs[2].imag], r3=[t.coeffs[3].real, t.coeffs[3].imag], G=G))
    print(f"q={q:>5} κ(1/q)={t.coeffs[2]:+.7f}  r3={t.coeffs[3]:+.6f}  G={G:.6f}", flush=True)
json.dump(rows, open("data/p1_limit.json", "w"))
