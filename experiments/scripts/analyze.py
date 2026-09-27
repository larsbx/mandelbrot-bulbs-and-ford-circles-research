"""Tables from data/taylor_q<q>.json (+ data/index_q<q>.json when present)."""
from paths import DATA
import json, sys, os
import numpy as np

def load(q):
    rows = json.load(open(DATA / f"taylor_q{q}.json"))
    for r in rows:
        r["r"] = np.array([complex(a, b) for a, b in r["r"]])
    idx = {}
    if os.path.exists(DATA / f"index_q{q}.json"):
        idx = {d["p"]: complex(d["re"], d["im"]) for d in json.load(open(DATA / f"index_q{q}.json"))}
    return sorted(rows, key=lambda r: r["xstar"]), idx

if __name__ == "__main__":
    for q in map(int, sys.argv[1:]):
        rows, idx = load(q)
        print(f"\n=== q={q}:  x*  p  | r2=κ (Re,Im) | r3 | r4 | G_ant  G_ant(series)  G_cen  G_cen(series) | tail | ι check")
        for r in rows:
            k = r["r"][2]
            chk = "" if r["p"] not in idx else f"{abs(k - (idx[r['p']] - 0.5) / q):.1e}"
            print(f"{r['xstar']:.4f} {r['p']:>4} | {k.real:+.5f} {k.imag:+.5f} | {r['r'][3].real:+.5f} {r['r'][3].imag:+.5f} | "
                  f"{r['r'][4].real:+.5f} {r['r'][4].imag:+.5f} | {r['G_ant']:.4f} {r['G_ant_series']:.4f}  {r['G_cen']:.4f} {r['G_cen_series']:.4f} | {r['tail']:.0e} | {chk}")
        e0 = max(abs(r["r"][0] - 1) for r in rows); e1 = max(abs(r["r"][1] + 1) for r in rows)
        print(f"max |r0-1| = {e0:.1e}, max |r1+1| = {e1:.1e}, max |G_ant - series| = {max(abs(r['G_ant']-r['G_ant_series']) for r in rows):.1e}")
