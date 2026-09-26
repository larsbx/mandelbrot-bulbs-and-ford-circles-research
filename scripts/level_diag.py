"""Diagonal test of the hierarchy on bulbs (paper Conjecture 6.3, V44): κ(p/q) for p/q = [0; N, …, N]
(k equal partial quotients), to be extrapolated in 1/N.  k = 2 is the control (limit κ̄₂ known, V42).
Uncertified normal form (Arb radii blow up for generic p): midpoints at 128 and 192 bits, the difference is the
reported error.  Writes data/level_diag.tsv: k, N, p, q, Re κ, Im κ, error."""
from multiprocessing import Pool
import mpmath as mp
from bulbford.cf import from_cf
from bulbford.normal_form_ball import normal_form_ball

JOBS = [(2, N) for N in range(16, 161, 8)] + [(3, N) for N in range(6, 41, 2)]


def run(job):
    k, N = job
    p, q = from_cf((0,) + (N,) * k)
    kap = normal_form_ball(p, q, prec=128, certified=False).kappa
    return k, N, p, q, kap.real.mid().str(25, radius=False), kap.imag.mid().str(25, radius=False), float(kap.rad())


if __name__ == "__main__":
    with Pool(3) as pool, open("data/level_diag.tsv", "w") as out:
        for row in pool.imap(run, sorted(JOBS, key=lambda j: -from_cf((0,) + (j[1],) * j[0])[1])):
            out.write("\t".join(map(str, row)) + "\n"); out.flush(); print(*row, flush=True)
