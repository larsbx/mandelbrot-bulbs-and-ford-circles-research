"""G and κ along canonical approaches x̃ = [0; tail, N, a] to a rational (fixed deep digit a, N → ∞):
is G linear in δ = |x̃ − p'/q'| along one family?  Families: ⅓⁻ [0;3,N,a], ⅓⁺ [0;2,1,N,a], ½⁻ [0;2,N,a], 0⁺ [0;N,a]."""
import json, sys
from bulbford.cf import from_cf, modinv
from bulbford.taylor import kappa_fft
FAM = {"1/3-": (lambda N, a: (0, 3, N, a), 1/3), "1/3+": (lambda N, a: (0, 2, 1, N, a), 1/3),
       "1/2-": (lambda N, a: (0, 2, N, a), 1/2), "0+": (lambda N, a: (0, N, a), 0.0)}
names = sys.argv[1:] or list(FAM)
rows = []
for name in names:
    f, x0 = FAM[name]
    for a in (16, 5):
        for N in (16, 24, 32, 48, 64, 96, 128, 192, 256, 384):
            pb, q = from_cf(f(N, a))
            if q > 20000: continue
            p = modinv(pb, q); t = kappa_fft(p, q); G = abs(t.solve(-1, 2.0)) / 2
            xt = (pb if pb <= q / 2 else pb - q) / q; d = abs(abs(xt) - x0)
            rows.append(dict(fam=name, a=a, N=N, p=p, q=q, xt=xt, delta=d, kappa=[t.coeffs[2].real, t.coeffs[2].imag], G=G))
            print(f"{name:<5} a={a:>2} N={N:>3} q={q:<6} δ={d:.2e} κ={t.coeffs[2]:+.5f} G={G:.5f}", flush=True)
json.dump(rows, open("data/canonical_profile.json", "w"))
