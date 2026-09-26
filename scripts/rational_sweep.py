"""Both one-sided constants κ^±_{p/q} = ι(𝒫^±_{p/q}) − ½ for all p/q with q ≤ QMAX, p ≤ q/2 (the others follow from
κ(1−x) = conj κ(x), which exchanges above and below).  Writes data/rational_sweep.tsv: p, q, Re/Im κ⁺, Re/Im κ⁻."""
import sys
from math import gcd
import mpmath as mp
from bulbford.horn_rational import horn_index_pq

QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
with open("data/rational_sweep.tsv", "w") as out:
    for q in range(1, QMAX + 1):
        for p in range(0 if q == 1 else 1, q // 2 + 1):
            if gcd(p, q) != 1:
                continue
            up, lo = (horn_index_pq(p, q, upper=u, dps=25) for u in (True, False))
            row = [p, q] + [mp.nstr(x, 16) for x in (mp.re(up), mp.im(up), mp.re(lo), mp.im(lo))]
            out.write("\t".join(map(str, row)) + "\n"); out.flush(); print(*row, flush=True)
