"""V51, V56: the root principle at a three-petal root.  Satellites of the period-3 component W_{1/3} (centre ≈ −0.1226 +
0.7449i, parametrised by the multiplier of its 3-cycle, dynamics.component) against the upper/lower parabolic
renormalizations of e^{2πi/3}z + z² at its root: angle [0;N,τ] (→ 0⁺, 𝒫⁺ with μ_τ) and 1 − [0;N,τ] (→ 1⁻, 𝒫⁻ with
μ̄_τ), [0;τ] = r/p.  V56 adds the tails with p_τ = 2, 3 on both sides, for κ and the bulb size G (no symmetry relates
the two sides of W_{1/3}, so the lower rows test the conjugate phase directly).  Writes data/period3_root.txt."""
import mpmath as mp
from bulbford.dynamics import component, bulb
from bulbford.taylor import taylor
from bulbford.extrapolate import power_fit
from bulbford.renorm import G_limits, bounded_p_limit
from bulbford.horn_rational import root, horn_index_pq

mp.mp.dps = 30
W3 = component("W1/3", 3, -0.1225611668766536 + 0.7448617666197442j)
NS_KAPPA, NS_G = (128, 256, 512, 1024, 2048), (256, 512, 1024, 2048, 4096)
# (quantity, tail (p, r) or None for τ = ∅, upper side?)
CASES = [("kappa", None, True), ("kappa", None, False), ("kappa", (2, 1), True),
         ("G", None, True), ("G", None, False),
         ("kappa", (2, 1), False), ("kappa", (3, 1), True), ("kappa", (3, 1), False),          # V56
         ("kappa", (3, 2), True), ("kappa", (3, 2), False),
         ("G", (2, 1), True), ("G", (2, 1), False), ("G", (3, 1), True), ("G", (3, 1), False)]


def angle(tail, upper, N):
    """(numerator, q) of [0;N,τ] (upper) or 1 − [0;N,τ] (lower); τ = ∅ gives 1/N."""
    p, r = tail or (1, 0)
    q = p * N + r
    return (p if upper else q - p), q


def label(qty, tail, upper):
    if tail is None:
        return f"{qty} {'1/q' if upper else '(q-1)/q'}"
    p, r = tail
    q = f"({p}N+{r})"
    return f"{qty} {p}/{q}" if upper else f"{qty} ({p}N-{p - r})/{q}"


def prediction(qty, tail, upper):
    germ = root(1, 3, upper)
    if qty == "kappa":
        return horn_index_pq(1, 3, upper=upper, dps=30) if tail is None else bounded_p_limit(*tail, dps=60, M=20, germ=germ)
    return G_limits(*(tail or (1, 1)), germ=germ)[0]


def bulbs(qty, tail, upper):
    if qty == "kappa":
        return [(q, complex(taylor(k, q, W3, r=1.5, N=64).coeffs[2])) for k, q in (angle(tail, upper, N) for N in NS_KAPPA)]
    return [(q, bulb(W3, k, q).G_ant) for k, q in (angle(tail, upper, N) for N in NS_G)]


if __name__ == "__main__":
    with open("data/period3_root.txt", "w") as out:
        for case in CASES:
            pred, rows = prediction(*case), bulbs(*case)
            f2, f3 = power_fit(rows, 2)[0], power_fit(rows, 3)[0]
            line = f"{label(*case)}\t{mp.nstr(pred, 12)}\t{mp.nstr(f3, 12)}\t{mp.nstr(abs(f3 - pred), 2)}\t{mp.nstr(abs(f3 - f2), 2)}"
            out.write(line + "\n"); out.flush(); print(line, flush=True)
