"""a(p/q) := [z^{q+1}] f_{λ₀}^q (leading parabolic coefficient) and b := [z^{2q+1}], vs x̃.  Ball arithmetic."""
import json, sys, cmath
from flint import acb, acb_series, arb, ctx
from bulbford.cf import coprime_numerators, modinv
from bulbford.index import _series_fq

def coeffs(p, q):
    old, oc = ctx.dps, ctx.cap
    try:
        ctx.dps = max(30, int(1.5 * q)); ctx.cap = 2 * q + 2
        s = _series_fq(p, q, 2 * q + 2).coeffs()
        loga = lambda z: (float(abs(z).log()), float(z.arg()))          # (log|z|, arg z) without overflow
        return loga(s[q + 1]), (loga(s[2 * q + 1]) if len(s) > 2 * q + 1 else (float('-inf'), 0.0))
    finally:
        ctx.dps, ctx.cap = old, oc

for q in map(int, sys.argv[1:] or ["59"]):
    rows = []
    for p in coprime_numerators(q):
        if p > q // 2: continue
        (la, aa), (lb, ab) = coeffs(p, q); pb = modinv(p, q); xt = (pb if pb <= q / 2 else pb - q) / q
        rows.append(dict(p=p, q=q, xt=xt, loga=la, arga=aa, logb=lb, argb=ab))
    rows.sort(key=lambda r: abs(r["xt"]))
    json.dump(rows, open(f"data/leading_coeff_q{q}.json", "w"))
    print(f"q={q}:   x̃      p    log|a|/q   arg(a)/2π    log|b|/q    |b/a²|/q   arg(b/a²)/2π")
    for r in rows:
        print(f"  {r['xt']:+.4f} {r['p']:>4}  {r['loga']/q:+.4f}   {r['arga']/(2*cmath.pi):+.4f}     {r['logb']/q:+.4f}    {cmath.exp(r['logb']-2*r['loga'])/q:8.4g}   {((r['argb']-2*r['arga'])/(2*cmath.pi)+0.5)%1-0.5:+.4f}")
