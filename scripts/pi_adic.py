"""V58: π-adic structure of the exact invariants at prime powers q = ℓ^v, π = 1 − ζ (the unique prime above ℓ in
ℚ(ζ_q), totally ramified, residue field 𝔽_ℓ).  For a = [z^{q+1}] f^q and X = q·ι·a² (bulbford/modular.py) the
valuation and leading residue (pi_adic); ι follows from ι = X/(q a²) with res(q/π^{vφ(q)}) = (−1)^v (ℓ odd; Wilson:
ℓ = ∏(1 − ζ^j), ∏ j ≡ −1).  For ℓ odd also v(R), R = ι − 1/((2 − ζ)^q − 1): the second fixed point z* = 1 − ζ of
ζz + z² (multiplier 2 − ζ) carries the leading π-adic term of ι.  Also every prime ℓ < 140.  Writes data/pi_adic.txt."""
import flint
from bulbford.modular import exact_invariants, pi_adic

PRIME_POWERS = [(3, 4), (5, 3), (7, 2), (11, 2), (13, 1), (2, 6)]
PRIMES = [l for l in range(3, 140) if flint.fmpz(l).is_prime()]


def row(l, v):
    q = l ** v
    e = exact_invariants(q)
    (va, ra), (vx, rx) = pi_adic(e.a, q), pi_adic(e.X, q)
    vi = vx - v * len(e.a) - 2 * va
    ri = rx * pow((-1) ** v % l if l > 2 else 1, -1, l) * pow(ra, -2, l) % l
    if l == 2:
        return q, va, ra, vx, rx, vi, ri, "-"
    phi, P = flint.fmpz_poly.cyclotomic(q), lambda c: flint.fmpz_poly([int(x) for x in c])
    m1 = (flint.fmpz_poly([2, -1]) ** q - 1) % phi                      # (2 − ζ)^q − 1
    Y = (P(e.X) * m1 - q * P(e.a) ** 2) % phi                             # q a² ((2 − ζ)^q − 1) R
    vr = pi_adic([int(c) for c in Y.coeffs()], q)[0] - v * len(e.a) - 2 * va - pi_adic([int(c) for c in m1.coeffs()], q)[0]
    return q, va, ra, vx, rx, vi, ri, vr


with open("data/pi_adic.txt", "w") as out:
    def emit(s):
        out.write(s + "\n"); out.flush(); print(s, flush=True)
    emit("q\tv(a)\tres(a)\tv(X)\tres(X)\tv(iota)\tres(iota·π^-v)\tv(R)\tlaw (odd l): v(a)=q/l, res(a)=-1, v(iota)=-(q+q/l), res=1; v(R)=l-1 (q=l), -2φ(q)/l (v≥2)")
    for l, vmax in PRIME_POWERS:
        for v in range(1, vmax + 1):
            emit("\t".join(map(str, row(l, v))))
    for l in PRIMES:
        emit("\t".join(map(str, row(l, 1))))
