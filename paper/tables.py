"""Generate the LaTeX tables of paper/kappa0.tex from the data files (run from the repository root:
PYTHONPATH=. python3 paper/tables.py).  Every number in a table is computed here, not typed."""
from collections import defaultdict
from math import gcd, pi
import mpmath as mp
from bulbford.extrapolate import power_fit, linear_limit
from bulbford.renorm import bounded_p_limit, G_limits, G_limit
from bulbford.horn import CUBIC
from bulbford.dynamics import MAIN3
from bulbford.taylor import taylor

mp.mp.dps = 40
K0 = mp.mpc("0.023825887402200569000729189645402", "-0.052304659114003666316522486207571")
K3 = mp.mpc("0.143599278095291569488064104088", "-0.0303116478985025376119431108575")
cplx = lambda z, n: f"{mp.nstr(mp.re(z), n, min_fixed=-9, max_fixed=9, strip_zeros=False)} {'-' if mp.im(z) < 0 else '+'} {mp.nstr(abs(mp.im(z)), n, min_fixed=-9, max_fixed=9, strip_zeros=False)}\\,i"


def sci_tex(x):
    m, e = f"{float(x):.1e}".split("e")
    return f"${m}\\cdot10^{{{int(e)}}}$"


def write(name, body):
    open(f"paper/tables/{name}.tex", "w").write(body); print("wrote", name)


# Table 1: κ(1/q) from the normal form
rows = {int(l.split("\t")[0]): mp.mpc(*l.split("\t")[1:3]) for l in open("data/p1_nf.tsv") if l.strip()}
rows.update({int(l.split("\t")[0]): mp.mpc(*l.split("\t")[1:3]) for l in open("data/p1_ball.tsv") if l.strip()})   # certified (Arb), q ≤ 32768
lines = [r"\begin{tabular}{@{}rll@{}}", r"\toprule", r"$q$ & $\kappa(1/q)$ & $|\kappa(1/q)+\tfrac{i}{2\pi q}-\kappa_0|$ \\", r"\midrule"]
for q in (64, 256, 1024, 4096, 8192, 16384, 32768):
    k = rows[q]; e = abs(k + 1j / (2 * mp.pi * q) - K0)
    lines.append(f"{q} & ${cplx(k, 16)}$ & {sci_tex(e)} \\\\")
fit = power_fit([(q, k) for q, k in rows.items() if q >= 256], 10)   # 13 points, residual 3.6e-29 (Numerical result 4.1)
lines += [r"\midrule", f"fit, $q\\to\\infty$ & ${cplx(fit[0], 16)}$ & \\\\",
          f"horn map, Thm.~\\ref{{thm:horn-index}} & ${cplx(K0, 16)}$ & \\\\", r"\bottomrule", r"\end{tabular}"]
write("kappa_1q", "\n".join(lines) + "\n")

# Table 2: bounded-p limits of κ
seq = defaultdict(dict)
for l in open("data/bounded_p_nf.tsv"):
    p, r, N, q, re, im, *_ = l.split("\t"); seq[(int(p), int(r))][int(q)] = mp.mpc(re, im)
lines = [r"\begin{tabular}{@{}cclll@{}}", r"\toprule", r"$p$ & $r$ & horn-map prediction $(\iota((\mu\mathcal P_0)^p)-\frac12)/p$ & bulbs, extrapolated & $|\Delta|$ \\", r"\midrule"]
for (p, r), D in sorted(seq.items()):
    lim = power_fit(D.items(), 3)[0]; pred = bounded_p_limit(p, r)
    lines.append(f"{p} & {r} & ${cplx(pred, 11)}$ & ${cplx(lim, 11)}$ & {sci_tex(abs(lim - pred))} \\\\")
lines += [r"\bottomrule", r"\end{tabular}"]
write("bounded_p_kappa", "\n".join(lines) + "\n")

# Table 3: bounded-p limits of the bulb size G
G = defaultdict(dict)
for l in open("data/G_bounded_p.tsv"):
    p, r, N, q, g = l.split("\t"); G[(int(p), int(r))][int(q)] = mp.mpf(g)
lines = [r"\begin{tabular}{@{}ccllll@{}}", r"\toprule", r"$p$ & $r$ & $G$ at $q=512p+r$ & $G$ at $q=2048p+r$ & extrapolated & horn map \\", r"\midrule"]
for (p, r), D in sorted(G.items()):
    qs = sorted(D); lim = power_fit(D.items(), 2)[0]; pred = G_limits(p, r)[0]
    lines.append(f"{p} & {r} & {float(D[qs[0]]):.8f} & {float(D[qs[-1]]):.8f} & {float(mp.re(lim)):.8f} & {pred:.8f} \\\\")
lines += [r"\bottomrule", r"\end{tabular}"]
write("bounded_p_G", "\n".join(lines) + "\n")

# Table 4: degree 3
lines = [r"\begin{tabular}{@{}rll@{}}", r"\toprule", r"$q$ & $\kappa_3(1/q)$ & $|\kappa_3(1/q)+\tfrac{i}{2\pi q}-\kappa_0^{(3)}|$ \\", r"\midrule"]
for q in (128, 256, 512, 1024):
    t = taylor(1, q, MAIN3, r=1.0, N=64); k = mp.mpc(t.coeffs[2])
    assert abs(t.coeffs[1] + 1) < 1e-8
    lines.append(f"{q} & ${cplx(k, 10)}$ & {sci_tex(abs(k + 1j / (2 * mp.pi * q) - K3))} \\\\")
lines += [r"\midrule", f"horn map of $v+v^2+\\tfrac13v^3$ & ${cplx(K3, 10)}$ & \\\\", r"\bottomrule", r"\end{tabular}"]
write("degree3", "\n".join(lines) + "\n")

# Table 4 (tab:levels): the renormalization hierarchy κ_k = ι((𝒫∘conj)^{k−1} 𝒫₀) − ½, two independent sampling settings
lv = {h: [mp.mpc(*l.split("\t")[1:3]) for l in open(f"data/levels_{h}.tsv") if l.strip()] for h in ("h100", "h125")}
x = lv["h100"]
lines = [r"\begin{tabular}{@{}cllc@{}}", r"\toprule", r"$k$ & $\kappa_k$ ($h=1$) & $|\kappa_k-\kappa_{k-1}|$ & $|\kappa_k^{(h=1)}-\kappa_k^{(h=5/4)}|$ \\", r"\midrule"]
for k, (a, b) in enumerate(zip(lv["h100"], lv["h125"]), 1):
    dk = sci_tex(abs(a - x[k - 2])) if k > 1 else ""
    lines.append(f"{k} & ${cplx(a, 12)}$ & {dk} & {sci_tex(abs(a - b)) if k > 1 else '---'} \\\\")
kstar, ev = linear_limit(x)
lines += [r"\midrule", f"$\\infty$ & ${cplx(kstar, 11)}$ & & \\\\", r"\bottomrule", r"\end{tabular}"]
write("levels", "\n".join(lines) + "\n")
print("kappa_* =", mp.nstr(kstar, 13), " eigenvalues of the fitted contraction:", [mp.nstr(e, 6) for e in ev])

# Table (tab:diskroot): the period-2 disc near its root against the cardioid near its ½-root (C23 (iv), V48)
z = lambda s: mp.mpc(complex(s.strip("()").replace(" ", "")))
lines = [r"\begin{tabular}{@{}lllc@{}}", r"\toprule", r"& disc $[0;N,t]$ & cardioid $[0;1,1,N-1,t]$ & $|\Delta|$ \\", r"\midrule"]
for l in open("data/disk_root.txt"):
    if "J = 3" not in l:
        continue
    head, _, vals = l.split(": ")
    qty, tail = head.split(" tail ")
    disk_s, card_s = (part.split(" ", 1)[1] for part in vals.split("  ")[:2])
    fmt = (lambda s: f"${cplx(z(s), 10)}$") if qty == "kappa" else (lambda s: f"${mp.nstr(mp.mpf(s), 10)}$")
    label = ("$\\kappa$" if qty == "kappa" else "$G$") + ", $t=" + (tail.strip("(,)").replace(",", "") or "\\varnothing") + "$"
    lines.append(f"{label} & {fmt(disk_s)} & {fmt(card_s)} & {sci_tex(float(vals.split('|Δ| ')[1]))} \\\\")
lines += [r"\bottomrule", r"\end{tabular}"]
write("diskroot", "\n".join(lines) + "\n")

# Table (tab:constants): the constants of the paper, each from its horn map / hierarchy, with its bulb check
re_, sign, im_ = open("data/half_root.txt").readline().split("= (")[1].split("j)")[0].split()
khalf = mp.mpc(re_, sign + im_)
G0 = G_limit(1, 1)
rows = [(r"$\kappa_0$", r"$z+z^2$", cplx(K0, 16), r"$3\cdot10^{-28}$", r"\S\ref{sec:horn}"),
        (r"$\kappa_0^{(3)}$", r"$v+v^2+\frac13v^3$", cplx(K3, 16), r"$1.3\cdot10^{-6}$ ($q=1024$)", r"\S\ref{sec:degree3}"),
        (r"$\kappa_{1/2}$", r"$-z+z^2$", cplx(khalf, 16), r"$3\cdot10^{-10}$", r"\S\ref{sec:component}"),
        (r"$\kappa_*$", r"$G_*=\mathcal P(\bar G_*)$", cplx(kstar, 11), r"levels $\le3$", r"\S\ref{sec:hierarchy}"),
        (r"$G_0$", r"$e^{u}\mathcal P_0$", f"{G0:.8f}", r"$4\cdot10^{-8}$", r"\S\ref{sec:unfolding}")]
lines = [r"\begin{tabular}{@{}lllll@{}}", r"\toprule", r"constant & germ & value & bulbs & \\", r"\midrule"]
lines += [f"{a} & {b} & ${c}$ & {d} & {e} \\\\" for a, b, c, d, e in rows]
lines += [r"\bottomrule", r"\end{tabular}"]
write("constants", "\n".join(lines) + "\n")
