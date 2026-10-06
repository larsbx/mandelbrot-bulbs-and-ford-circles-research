"""V57: exact a and X = q·ι·a² in ℚ(ζ_q) by multimodular arithmetic (bulbford/modular.py) for every q in 2..QMAX
and the extra q given.  Per q: φ(q), primes used, seconds, coordinate bits of a and X, denom(X) (conjecture C22: 1),
d(q) = denom(ι·a²) = denom(X/q), and the relative error of σ₁(X) against the floating normal form.  Writes
data/modular_invariants.tsv; the exact coordinates for q ≤ 32 go to data/modular_exact.txt."""
import sys, time
from fractions import Fraction
from math import lcm
import mpmath as mp
from bulbford.modular import exact_invariants, embed
from bulbford.normal_form import normal_form

QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
QS = list(range(2, QMAX + 1)) + [int(a) for a in sys.argv[2:]]

with open("data/modular_invariants.tsv", "w") as out, open("data/modular_exact.txt", "w") as ex:
    out.write("q\tphi\tprimes\tseconds\tbits_a\tbits_X\tdenom_X\tdenom_iota_a2\trel_err\n")
    for q in QS:
        t = time.time()
        e = exact_invariants(q)
        dt = time.time() - t
        nf = normal_form(1, q)
        with mp.workdps(q + 40):
            ref = q * nf.iota * nf.a ** 2
            err = abs(embed(e.X, q) - ref) / abs(ref)
        row = (q, len(e.a), e.primes, f"{dt:.2f}", max(abs(x) for x in e.a).bit_length(),
               max(abs(x.numerator) for x in e.X).bit_length(), lcm(*(x.denominator for x in e.X)),
               lcm(*(Fraction(x, q).denominator for x in e.X)), mp.nstr(err, 2))
        out.write("\t".join(map(str, row)) + "\n"); out.flush(); print(*row, flush=True)
        if q <= 32:
            ex.write(f"q={q}\ta={list(e.a)}\tX={[str(x) for x in e.X]}\n"); ex.flush()
