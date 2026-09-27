"""Rigorous ι_{p/q} (ball arithmetic) for p ≤ q/2; writes data/index_q<q>.json with mid/rad."""
from paths import DATA
import json, sys, time
from bulbford.cf import coprime_numerators, xstar
from bulbford.index import index

if __name__ == "__main__":
    for q in map(int, sys.argv[1:]):
        t0 = time.time(); rows = []
        for p in coprime_numerators(q):
            if p > q // 2: continue
            r = index(p, q)
            rows.append(dict(p=p, q=q, xstar=float(xstar(p, q)), re=float(r.real.mid()), im=float(r.imag.mid()), rad=float(r.rad())))
        json.dump(rows, open(DATA / f"index_q{q}.json", "w"))
        print(f"q={q}: {len(rows)} values, max rad {max(r['rad'] for r in rows):.1e}, {time.time()-t0:.1f}s", flush=True)
