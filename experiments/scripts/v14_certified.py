"""V14 with both proxies certified: Δ = ||u_a|/2 − G_ant| from the P10 antipode boxes.

Reads experiments/data/antipode_certificates_v14.json (replayed by the tests), encloses |u_a|/2 on each certified
c box with Arb (bulbford.antipode.u_half_abs), and writes experiments/data/v14_certified.json: per-p brackets for
G_ant, |u_a|/2 and Δ, and per-q brackets for max_p Δ and q²·max_p Δ. `--check` exits 1 on drift.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

from bulbford.antipode import u_half_abs
from bulbford.certify import I, box_from_numerators
from paths import DATA

SOURCE = DATA / "antipode_certificates_v14.json"
OUT = DATA / "v14_certified.json"


def outward(x: Fraction, digits: int, up: bool) -> str:
    """x as a decimal rounded away from the enclosure's inside (down for a lower end, up for an upper)."""
    scaled = x * 10**digits
    k = -((-scaled.numerator) // scaled.denominator) if up else scaled.numerator // scaled.denominator
    sign, k = ("-", -k) if k < 0 else ("", k)
    return f"{sign}{k // 10**digits}.{k % 10**digits:0{digits}d}"


def abs_interval(d: I) -> I:
    if d.lo > 0:
        return d
    if d.hi < 0:
        return I(-d.hi, -d.lo)
    return I(Fraction(0), max(-d.lo, d.hi))


def row(r: dict) -> dict:
    g = I(*(Fraction(x) for x in r["G_ant_bounds"]))
    u = u_half_abs(r["p"], r["q"], box_from_numerators(r["c_box_numerators"], r["prec"]), r["prec"])
    d = abs_interval(u - g)
    return {"p": r["p"], "q": r["q"], "G_ant": [str(g.lo), str(g.hi)], "u_a_half": [str(u.lo), str(u.hi)],
            "delta": [str(d.lo), str(d.hi)]}


def summary(rows: list[dict], q: int) -> dict:
    ds = [I(*(Fraction(x) for x in r["delta"])) for r in rows if r["q"] == q]
    top = I(max(d.lo for d in ds), max(d.hi for d in ds))
    arg = max((r for r in rows if r["q"] == q), key=lambda r: Fraction(r["delta"][1]))["p"]
    return {"q": q, "max_delta": [str(top.lo), str(top.hi)], "argmax_p": arg,
            "q2_max_delta": [str(top.lo * q * q), str(top.hi * q * q)]}


def render() -> str:
    records = json.loads(SOURCE.read_text())["certificates"]
    rows = [row(r) for r in records]
    qs = sorted({r["q"] for r in rows})
    doc = {"schema": "bulbford-v14-certified/v1", "source": SOURCE.name,
           "definition": "delta = | |u_a|/2 - G_ant |, u_a = q^2 log(lambda_a/lambda_0), lambda_a = 1 - sqrt(1 - 4 c_ant)",
           "summary": [summary(rows, q) for q in qs], "rows": rows}
    return json.dumps(doc, indent=1) + "\n"


if __name__ == "__main__":
    text = render()
    if "--check" in sys.argv[1:]:
        sys.exit(0 if OUT.read_text() == text else 1)
    OUT.write_text(text)
    for s in json.loads(text)["summary"]:
        lo, hi = (Fraction(x) for x in s["q2_max_delta"])
        print(f"q={s['q']}: max Δ at p={s['argmax_p']}, q²Δ ∈ [{outward(lo, 12, up=False)}, {outward(hi, 12, up=True)}]")
