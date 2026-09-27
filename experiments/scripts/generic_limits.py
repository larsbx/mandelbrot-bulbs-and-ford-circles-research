"""Generic (p → ∞) one-sided limits: x̃ = tail + [a1] with a1 large, N growing.  Families:
⅓⁺: [0;2,1,N,a1];  ⅓⁻: [0;3,N,a1];  0⁺: [0;N,a1];  ½⁻: [0;2,N,a1]."""
from paths import DATA
import json, sys
from bulbford.cf import from_cf, modinv
from bulbford.taylor import kappa_fft
FAM = {"1/3+": lambda N, a: (0, 2, 1, N, a), "1/3-": lambda N, a: (0, 3, N, a),
       "0+": lambda N, a: (0, N, a), "1/2-": lambda N, a: (0, 2, N, a)}
names = sys.argv[1:] or list(FAM)
rows = []
for name in names:
    for N in (64, 128, 256):
        for a in (16, 64):
            pb, q = from_cf(FAM[name](N, a))
            if q > 20000: continue
            p = modinv(pb, q); t = kappa_fft(p, q); G = abs(t.solve(-1, 2.0)) / 2
            xt = (pb if pb <= q / 2 else pb - q) / q
            rows.append(dict(fam=name, N=N, a1=a, p=p, q=q, xt=xt, kappa=[t.coeffs[2].real, t.coeffs[2].imag], G=G))
            print(f"{name:<5} N={N:>3} a1={a:>2} p/q={p}/{q:<6} x̃={xt:+.5f} κ={t.coeffs[2]:+.5f} G={G:.5f}", flush=True)
json.dump(rows, open(DATA / "generic_limits.json", "w"))
