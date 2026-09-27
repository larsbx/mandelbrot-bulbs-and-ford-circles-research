import numpy as np, math
from math import gcd, pi, sin, log

def iterate(c, z, q):
    """return f_c^q(z), d/dz, d/dc"""
    A, B = 1.0+0j, 0j
    for _ in range(q):
        A, B = 2*z*A, 2*z*B + 1
        z = z*z + c
    return z, A, B

def center(c0, q):
    c = c0
    for _ in range(60):
        z, A, B = iterate(c, 0j, q)
        step = z / B
        c -= step
        if abs(step) < 1e-15: break
    return c

def cycle_point(c, z0, q):
    z = z0
    for _ in range(60):
        w, A, B = iterate(c, z, q)
        step = (w - z)/(A - 1)
        z -= step
        if abs(step) < 1e-15: break
    return z

def multiplier(c, z, q):
    return iterate(c, z, q)[1]

def antipode(c_root, c_cen, q):
    """boundary point with multiplier -1, continued from center along the root->center ray"""
    d = c_cen - c_root
    z = 0j
    t = 0.0
    c = c_cen
    for k in range(1, 25):           # continuation to t ~ 1.2 (predicted antipode at t=1)
        t = 1.2*k/24
        c = c_root + (1+t)*d
        z = cycle_point(c, z, q)
    # Newton in c on rho(c) + 1 = 0, cycle tracked by Newton in z
    for _ in range(60):
        z = cycle_point(c, z, q)
        r = multiplier(c, z, q)
        h = 1e-7*abs(d)
        z2 = cycle_point(c+h, z, q); r2 = multiplier(c+h, z2, q)
        z3 = cycle_point(c-h, z, q); r3 = multiplier(c-h, z3, q)
        dr = (r2 - r3)/(2*h)
        step = (r + 1)/dr
        c -= step
        if abs(step) < 1e-14*abs(d): break
    z = cycle_point(c, z, q)
    return c, multiplier(c, z, q)

if __name__=="__main__":
  pass
Q = 60
rows = []
for q in range(2, Q+1):
    for p in range(1, q):
        if gcd(p, q) != 1: continue
        lam0 = np.exp(2j*pi*p/q)
        c_root = lam0/2 - lam0**2/4
        guess = c_root + (1-lam0)*lam0/(2*q*q)
        c_cen = center(guess, q)
        c_ant, rho = antipode(c_root, c_cen, q)
        if abs(rho+1) > 1e-8: print("warn", p, q, rho)
        d_cen = 2*abs(c_cen - c_root)
        d_ant = abs(c_ant - c_root)
        pred = 2*sin(pi*p/q)/q**2
        rows.append((p, q, d_cen, d_ant, pred))

print(f"{'p/q':>7} {'2|cen-root|':>12} {'|ant-root|':>12} {'pred':>12} {'ratio_cen':>10} {'ratio_ant':>10}")
for p,q,dc,da,pr in rows:
    if q <= 7 or (q in (20,40,60) and p in (1, q//3, q//2 if q%2 else q//2+1)):
        print(f"{p:>3}/{q:<3} {dc:12.6f} {da:12.6f} {pr:12.6f} {dc/pr:10.5f} {da/pr:10.5f}")

print("\nper-q: mean ratio (antipode), mean of q*(ratio-1) [c1 estimate], spread over p")
for q in (5,10,20,30,40,50,60):
    rs = [(p/q, da/pr) for p,qq,dc,da,pr in rows if qq==q]
    r = np.array([x[1] for x in rs])
    print(f"q={q:<3} mean ratio={r.mean():.5f}  q*(ratio-1): mean={q*(r-1).mean():.4f} min={q*(r-1).min():.4f} max={q*(r-1).max():.4f}")

print("\nc1(t) estimate at q=60, sampled t=p/q:")
for p,q,dc,da,pr in rows:
    if q==60 and p in (1,7,13,19,23,29,31,37,41,47,53,59):
        print(f"  t={p/q:.4f}  q*(ratio-1)={q*(da/pr-1):+.4f}")

print("\nSum of diameters vs (24/pi^3) log Q:")
K = 24/pi**3
S = 0.0; prevS=None; prevL=None
for q in range(2, Q+1):
    S += sum(da for p,qq,dc,da,pr in rows if qq==q)
    if q in (10,20,30,40,50,60):
        L = log(q)
        slope = (S-prevS)/(L-prevL) if prevS is not None else float('nan')
        print(f"Q={q:<3} sum d={S:.5f}  K log Q={K*L:.5f}  local slope dS/dlogQ={slope:.4f}  (K={K:.4f})")
        prevS, prevL = S, L
