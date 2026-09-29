"""Dense κ(p/q) via FFT for all p ≤ q/2 (double precision, N=64). Writes experiments/data/kappa_q<q>.json.

    PYTHONPATH=kernel python3 experiments/scripts/kappa_sweep.py [--jobs J] q [q ...]
"""
from paths import DATA
import argparse, json, os, time
from multiprocessing import Pool
from bulbford.cf import coprime_numerators, modinv
from bulbford.taylor import kappa_fft


def row(pq: tuple[int, int]) -> dict:
    p, q = pq
    t = kappa_fft(p, q); pb = modinv(p, q); ua = t.solve(-1, 2.0)
    return dict(p=p, q=q, xt=(pb if pb <= q / 2 else pb - q) / q, kappa=[t.coeffs[2].real, t.coeffs[2].imag],
                r3=[t.coeffs[3].real, t.coeffs[3].imag], G=abs(ua) / 2)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("qs", type=int, nargs="+")
    ap.add_argument("--jobs", type=int, default=1)
    args = ap.parse_args()
    for q in args.qs:
        t0 = time.time()
        output = DATA / f"kappa_q{q}.json"
        checkpoint = DATA / f".kappa_q{q}.jsonl"
        done = {}
        if checkpoint.exists():
            for line in checkpoint.read_text().splitlines():
                r = json.loads(line); done[r["p"]] = r
        work = [(p, q) for p in coprime_numerators(q) if p <= q // 2 and p not in done]
        with Pool(args.jobs) as pool:
            with checkpoint.open("a") as cp:
                for r in pool.imap_unordered(row, work, chunksize=1):
                    done[r["p"]] = r
                    cp.write(json.dumps(r) + "\n"); cp.flush()
        rows = [done[p] for p in sorted(done)]
        with output.open("w") as f:
            json.dump(rows, f)
        os.unlink(checkpoint)
        print(f"q={q}: {len(rows)} rows, {time.time()-t0:.0f}s, jobs={args.jobs}", flush=True)
