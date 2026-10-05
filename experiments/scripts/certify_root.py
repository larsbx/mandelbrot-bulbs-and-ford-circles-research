"""P14: the certified centre's component has the satellite root c_root(p/q) on its boundary.

For each bulb: an untrusted float survey picks the disk D around z₀ = λ₀/2 (radius r below the nearest other
periodic point of f^q at the root), D′ (r′ + r′² < r), and c₁ on the chord from c_root towards the P3 centre
close enough to the root that the cycle born there lies in D′. Then bulbford.root certifies the segment
c_root → c₁ and bulbford.continuation certifies the path centre → c₁ ending inside D′. Writes
experiments/data/root_certificates.json; `--check` recomputes everything and exits 1 on drift.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from multiprocessing import Pool

import numpy as np
from numpy.polynomial import polynomial as npoly
from flint import acb, ctx

from bulbford.antipode import lambda_box, root_box
from bulbford.certify import box_from_numerators
from bulbford.continuation import _acb, _arb_exact, continue_centre_to_point
from bulbford.index import lambda_ball
from bulbford.root import WORK_PREC, certify_root_side

from paths import DATA

OUT = DATA / "root_certificates.json"
Q_MAX = 4  # the recertified corpus; the previous angle-based integration had only partial q = 5 coverage


def _dyadic(x: float, bits: int = 40) -> Fraction:
    return Fraction(round(x * 2**bits), 2**bits)


def survey(p: int, q: int, centre: complex) -> dict:
    """Untrusted choice of r, r′ and τ (c₁ = c_root + τ(centre − c_root)); certificates decide."""
    lam = complex(lambda_ball(p, q, WORK_PREC))  # untrusted midpoint; certified checks decide acceptance
    z0, croot = lam / 2, lam / 2 - lam * lam / 4
    w = np.array([0, 1], dtype=complex)
    for _ in range(q):
        w = npoly.polyadd(npoly.polymul(w, w), [croot])
    dist = sorted(abs(npoly.polyroots(npoly.polysub(w, [0, 1])) - z0))
    r = 0.8 * dist[q + 1]  # the (q+2)-th root: past the (q+1)-fold cluster at z₀
    r_in = min((-1 + (1 + 4 * 0.9 * r) ** 0.5) / 2, 0.45)  # r′ + r′² = 0.9 r, and 2r′ < 1 for the log
    tau = 0.5
    while tau > 2**-24:
        c1 = croot + tau * (centre - croot)
        z = 0j  # the superattracting point at the centre, followed down the chord to c₁
        for k in range(1, 801):
            c = centre + (c1 - centre) * k / 800
            for _ in range(40):
                f, a = z, 1
                for _ in range(q):
                    a, f = 2 * f * a, f * f + c
                z -= (f - z) / (a - 1)
        orbit = [z]
        for _ in range(q - 1):
            orbit.append(orbit[-1] ** 2 + c1)
        if max(abs(o - z0) for o in orbit) < 0.8 * r_in:
            break
        tau /= 2
    return {"m": z0, "r": r, "r_in": r_in, "tau": tau}


def record(args: tuple) -> dict:
    old = ctx.prec
    try:
        ctx.prec = WORK_PREC
        return _record(args)
    finally:
        ctx.prec = old


def _record(args: tuple) -> dict:
    p, q, centre = args
    cmid = complex(float(centre.mid[0]), float(centre.mid[1]))
    s = survey(p, q, cmid)
    c_root = _acb(root_box(lambda_box(p, q)))
    # c₁ exact dyadic on the chord from mid(c_root) towards the centre
    crm = complex(float(c_root.real.mid()), float(c_root.imag.mid()))
    c1f = crm + s["tau"] * (cmid - crm)
    c1 = (_dyadic(c1f.real, 60), _dyadic(c1f.imag, 60))
    m = (_dyadic(s["m"].real), _dyadic(s["m"].imag))
    r, r_in = _dyadic(s["r"]), _dyadic(s["r_in"])
    return replay({"p": p, "q": q, "disk_centre": [str(m[0]), str(m[1])],
                   "r": str(r), "r_inner": str(r_in), "c1": [str(c1[0]), str(c1[1])]}, centre)


def replay(row: dict, centre) -> dict:
    """Recheck the stored exact disk and endpoint without consulting the untrusted survey."""
    old = ctx.prec
    try:
        ctx.prec = WORK_PREC
        p, q = row["p"], row["q"]
        m = tuple(Fraction(x) for x in row["disk_centre"])
        c1 = tuple(Fraction(x) for x in row["c1"])
        r, r_in = Fraction(row["r"]), Fraction(row["r_inner"])
        m_b, c1_b = acb(_arb_exact(m[0]), _arb_exact(m[1])), acb(_arb_exact(c1[0]), _arb_exact(c1[1]))
        c_root = _acb(root_box(lambda_box(p, q)))
        side = certify_root_side(p, q, c_root, c1_b, m_b, _arb_exact(r), _arb_exact(r_in))
        path = continue_centre_to_point(p, q, centre, c1, (m_b, _arb_exact(r_in)))
        return {"p": p, "q": q, "accepted": side.accepted and path.accepted, "disk_centre": list(row["disk_centre"]),
                "r": row["r"], "r_inner": row["r_inner"], "c1": list(row["c1"]), "root_side_boxes": side.boxes,
                "root_side": side.reason, "path_pieces": path.pieces, "path": path.reason}
    finally:
        ctx.prec = old


def jobs(q_max: int = Q_MAX) -> list[tuple]:
    cen = json.loads((DATA / "center_certificates.json").read_text())["certificates"]
    return [(r["p"], r["q"], box_from_numerators(r["box_numerators"], r["prec"]))
            for r in sorted(cen, key=lambda r: (r["q"], r["p"])) if r["q"] <= q_max]


def render(q_max: int = Q_MAX, workers: int = 4) -> str:
    with Pool(workers) as pool:
        rows = pool.map(record, jobs(q_max), chunksize=1)
    doc = {"schema": "bulbford-root-certificates/v1",
           "claim": "the hyperbolic component of the certified centre has c_root(p/q) on its boundary",
           "arithmetic": "Arb balls and algebraic-root contour means with certified Laurent error bounds",
           "certificates": rows}
    return json.dumps(doc, indent=1) + "\n"


if __name__ == "__main__":
    q_max = int(next((a.split("=")[1] for a in sys.argv[1:] if a.startswith("--q-max=")), Q_MAX))
    text = render(q_max)
    if "--check" in sys.argv[1:]:
        sys.exit(0 if OUT.read_text() == text else 1)
    OUT.write_text(text)
    rows = json.loads(text)["certificates"]
    print(f"wrote {OUT.name}: {sum(r['accepted'] for r in rows)}/{len(rows)} accepted")
    for r in rows:
        if not r["accepted"]:
            print("  not accepted:", r["p"], r["q"], r["root_side"], "|", r["path"])
