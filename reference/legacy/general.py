import numpy as np
from math import pi, sin, gcd, sqrt

def it(c, z, q, d):
    A, B = 1+0j, 0j
    for _ in range(q):
        w = d*z**(d-1)
        A, B = w*A, w*B + 1
        z = z**d + c
    return z, A, B

def center(c, q, d):
    for _ in range(80):
        z, A, B = it(c, 0j, q, d); s = z/B; c -= s
        if abs(s) < 1e-15: break
    return c

def cyc(c, z, q, d):
    for _ in range(80):
        w, A, B = it(c, z, q, d); s = (w-z)/(A-1); z -= s
        if abs(s) < 1e-15: break
    return z

def antipode(c_root, c_cen, q, d):
    dd = c_cen - c_root; z = 0j; c = c_cen
    for k in range(1, 25):
        c = c_root + (1 + 1.2*k/24)*dd; z = cyc(c, z, q, d)
    for _ in range(80):
        z = cyc(c, z, q, d); r = it(c, z, q, d)[1]; h = 1e-7*abs(dd)
        r2 = it(c+h, cyc(c+h, z, q, d), q, d)[1]; r3 = it(c-h, cyc(c-h, z, q, d), q, d)[1]
        s = (r+1)/((r2-r3)/(2*h)); c -= s
        if abs(s) < 1e-14*abs(dd): break
    return c

def G(family, p, q):
    lam0 = np.exp(2j*pi*p/q); eps = 1/q**2
    if family == 'main2':
        cl = lambda l: l/2 - l**2/4; d = 2; k = 1
        pred = 2*sin(pi*p/q)/q**2
    elif family == 'disk2':                      # period-2 disk, multiplier 4(c+1)
        cl = lambda l: l/4 - 1; d = 2; k = 2
        pred = 2*(1/4)/q**2
    elif family == 'main3':                      # z^3+c main component, param by z=sqrt(l/3)
        cl = lambda l: (lambda z: z - z**3)(np.sqrt(l/3)); d = 3; k = 1
        pred = 2*abs(1-lam0)/(6*abs(np.sqrt(lam0/3)))/q**2   # 2|c'(lam0)|/q^2
    c_root = cl(lam0); c_cen = center(cl(lam0*(1+eps)), k*q, d)
    c_ant = antipode(c_root, c_cen, k*q, d)
    return abs(c_ant - c_root)/pred, 2*abs(c_cen-c_root)/pred

if __name__ == "__main__":
    q = 59
    rows = []
    for p in range(1, q//2+1):
        pb = pow(p, -1, q); xs = min(pb, q-pb)/q
        g2, _ = G('main2', p, q); gd, _ = G('disk2', p, q); g3, _ = G('main3', p, q)
        rows.append((xs, p, g2, gd, g3))
    rows.sort()
    print("q=59, sorted by x*=||p^-1/q||.   x*     p    G_main2  G_disk2  G_main3")
    for xs,p,g2,gd,g3 in rows:
        print(f"   {xs:.4f}  {p:>3}   {g2:.4f}   {gd:.4f}   {g3:.4f}")
    a = np.array([[r[2],r[3],r[4]] for r in rows])
    print("\ncorr(main2,disk2) =", np.corrcoef(a[:,0],a[:,1])[0,1].round(4))
    print("corr(main2,main3) =", np.corrcoef(a[:,0],a[:,2])[0,1].round(4))
    print("mean ratio disk2/main2 =", (a[:,1]/a[:,0]).mean().round(4), " sd", (a[:,1]/a[:,0]).std().round(4))
    print("mean ratio main3/main2 =", (a[:,2]/a[:,0]).mean().round(4), " sd", (a[:,2]/a[:,0]).std().round(4))
    # continuity: max jump between x*-adjacent entries vs total range
    g = a[:,0]; print("\nmain2: total range", (g.max()-g.min()).round(4), " max adjacent-in-x* jump", np.abs(np.diff(g)).max().round(4))
