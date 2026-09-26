"""Figures for paper/kappa0.tex.  Run from the repository root: PYTHONPATH=. python3 paper/figures.py
Palette: categorical slots 1–3 of the reference palette (validated: CVD ΔE 9.2, normal 27.6; aqua < 3:1 on white,
so every series is also direct-labelled and has its own marker shape)."""
import json, os
from math import gcd
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from bulbford.renorm import bounded_p_limit

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dd"
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9, "axes.edgecolor": MUTED,
                     "axes.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": INK,
                     "xtick.labelcolor": INK, "ytick.labelcolor": INK, "legend.frameon": False})
mp.mp.dps = 40
K0 = mp.mpc("0.023825887402200569000729189645402", "-0.052304659114003666316522486207571")
K3 = mp.mpc("0.1435992780952920", "-0.0303116478985025")


def style(ax):
    ax.grid(True, color=GRID, linewidth=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)


# Figure 1: convergence of κ(1/q) to κ₀ (degree 2) and κ₃(1/q) to κ₀⁽³⁾ (degree 3)
merged = {int(r[0]): mp.mpc(r[1], r[2]) for r in (l.split("\t") for l in open("data/p1_nf.tsv") if l.strip())}
merged.update({int(r[0]): mp.mpc(r[1], r[2]) for r in (l.split("\t") for l in open("data/p1_ball.tsv") if l.strip())})
q = np.array(sorted(merged)); k = [merged[qq] for qq in q]
raw = np.array([float(abs(x - K0)) for x in k])
cor = np.array([float(abs(x + 1j / (2 * mp.pi * qq) - K0)) for x, qq in zip(k, q)])
q3 = np.array([128, 256, 512, 1024])
e3 = np.array([8.5e-5, 2.1e-5, 5.3e-6, 1.3e-6])      # data/main3_p1.txt (|κ₃ + i/(2πq) − κ₀⁽³⁾|)
fig, ax = plt.subplots(figsize=(4.8, 3.2)); style(ax)
ax.loglog(q, raw, "-o", color=BLUE, lw=2, ms=4.5, mec="white", mew=1.0, label=r"$|\kappa(1/q)-\kappa_0|$")
ax.loglog(q, cor, "-s", color=ORANGE, lw=2, ms=4.5, mec="white", mew=1.0, label=r"$|\kappa(1/q)+\frac{i}{2\pi q}-\kappa_0|$")
ax.loglog(q3, e3, "-^", color=AQUA, lw=2, ms=5, mec="white", mew=1.0, label=r"degree 3: $|\kappa_3(1/q)+\frac{i}{2\pi q}-\kappa_0^{(3)}|$")
for x0, y0, sl, lab in ((64, raw[0] * 1.9, -1, r"$\propto q^{-1}$"), (64, cor[0] * 0.28, -2, r"$\propto q^{-2}$")):
    xs = np.array([x0, 32768]); ax.loglog(xs, y0 * (xs / x0) ** sl, color=MUTED, lw=0.6)
    xl = 900 if sl == -1 else 3000
    ax.text(xl, y0 * (xl / x0) ** sl * (1.6 if sl == -1 else 0.45), lab, color=MUTED, fontsize=8)
ax.text(q[-1] * 1.12, raw[-1], "raw", color=INK, va="center", fontsize=8)
ax.text(q[-1] * 1.12, cor[-1], "corrected", color=INK, va="center", fontsize=8)
ax.text(q3[-1] * 1.12, e3[-1], "degree 3", color=INK, va="center", fontsize=8)
ax.set_xlim(50, 1.3e5); ax.set_xlabel(r"$q$"); ax.set_ylabel("error")
ax.legend(loc="lower left", fontsize=7.5)
fig.tight_layout(); fig.savefig("paper/figures/convergence.pdf")
if os.environ.get("PNG_DIR"): fig.savefig(os.path.join(os.environ["PNG_DIR"], "convergence.png"), dpi=150)
plt.close(fig)

# Figure 2: all κ(p/q) at q = 2003 and the horn-map predictions
d = json.load(open("data/kappa_q2003.json"))
kk = np.array([complex(*r["kappa"]) for r in d]); kk = np.r_[kk, np.conj(kk)]
preds = [(1, 1)] + [(p, r) for p in range(2, 6) for r in range(1, p) if gcd(p, r) == 1]
pv = np.array([complex(bounded_p_limit(p, r)) for p, r in preds]); pv = np.r_[pv, np.conj(pv)]
fig, ax = plt.subplots(figsize=(3.6, 4.6)); style(ax)
ax.scatter(kk.real, kk.imag, s=4, color=BLUE, alpha=0.55, linewidths=0, label=r"$\kappa(p/q)$, $q=2003$, all $p$")
ax.scatter(pv.real, pv.imag, s=34, marker="D", color=ORANGE, edgecolors="white", linewidths=1.0,
           label=r"horn-map limits, bounded $p\leq 5$", zorder=3)
k0 = complex(K0)
ax.scatter([k0.real, k0.real], [k0.imag, -k0.imag], s=120, marker="*", color=AQUA, edgecolors="white", linewidths=1.0,
           label=r"$\kappa_0,\ \overline{\kappa_0}$ ($p=1$)", zorder=4)
ax.annotate(r"$\kappa_0$", (k0.real, k0.imag), xytext=(k0.real + 0.004, k0.imag - 0.001), color=INK, fontsize=9)
ax.annotate(r"$\overline{\kappa_0}$", (k0.real, -k0.imag), xytext=(k0.real + 0.004, -k0.imag - 0.001), color=INK, fontsize=9)
ax.set_xlabel(r"$\mathrm{Re}\,\kappa$"); ax.set_ylabel(r"$\mathrm{Im}\,\kappa$"); ax.set_aspect("equal")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), fontsize=7.5, markerscale=1.4, ncol=1)
fig.tight_layout(); fig.savefig("paper/figures/kappa_cloud.pdf")
if os.environ.get("PNG_DIR"): fig.savefig(os.path.join(os.environ["PNG_DIR"], "kappa_cloud.png"), dpi=150)
plt.close(fig)


def save(fig, name):
    fig.tight_layout(); fig.savefig(f"paper/figures/{name}.pdf")
    if os.environ.get("PNG_DIR"): fig.savefig(os.path.join(os.environ["PNG_DIR"], f"{name}.png"), dpi=150)
    plt.close(fig)


# Figure 3: the first-order correction is kinematic (data/kinematic.txt, degree 2, p = 1)
block, rows = None, []
for l in open("data/kinematic.txt"):
    if l.startswith("degree"):
        block = l.split(":")[0]
    elif block == "degree 2, p = 1, r = 1":
        head, tail = l.split("[")
        rows.append((int(head.split("=")[1].split(":")[0]), [float(x) for x in head.split(":")[1].split()],
                     [float(x) for x in tail.split("]")[0].split()]))
qk = np.array([r[0] for r in rows]); moeb = np.array([r[1] for r in rows]); plain = np.array([r[2] for r in rows])
fig, ax = plt.subplots(figsize=(4.8, 3.2)); style(ax)
for j, (col, mk, k) in enumerate(((BLUE, "o", 2), (ORANGE, "s", 3))):
    ax.loglog(qk, plain[:, j], "--", marker=mk, color=col, lw=1.4, ms=4.5, mfc="white", mew=1.2,
              label=rf"$|r_{k}-\rho_{k}|$")
    ax.loglog(qk, moeb[:, j], "-", marker=mk, color=col, lw=2, ms=4.5, mec="white", mew=1.0,
              label=rf"$|r_{k}-[u^{k}]\rho(\tilde u)|$")
for y0, sl, lab in ((plain[0, 0] * 2.2, -1, r"$\propto q^{-1}$"), (moeb[0, 1] * 0.3, -2, r"$\propto q^{-2}$")):
    xs = np.array([qk[0], qk[-1]]); ax.loglog(xs, y0 * (xs / qk[0]) ** sl, color=MUTED, lw=0.6)
    ax.text(qk[-1] * 1.1, y0 * (qk[-1] / qk[0]) ** sl, lab, color=MUTED, fontsize=8, va="center")
ax.set_xlim(50, 2400); ax.set_xlabel(r"$q$"); ax.set_ylabel(r"residual ($p=1$)")
ax.legend(loc="lower left", fontsize=7.5, ncol=2)
save(fig, "kinematic")

# Figure 4: the two hierarchies (conjugated: bulbs; holomorphic: Lanford–Yampolsky)
seq = lambda f: [mp.mpc(*l.split("\t")[1:3]) for l in open(f) if l.strip()]
fig, ax = plt.subplots(figsize=(4.8, 3.0)); style(ax)
for f, col, mk, lab in (("data/levels_h100.tsv", BLUE, "o", r"$G\mapsto\mathcal{P}(\bar G)$ (bulbs)"),
                        ("data/levels_holo100.tsv", ORANGE, "s", r"$G\mapsto\mathcal{P}(G)$")):
    x = seq(f); dk = [float(abs(b - a)) for a, b in zip(x, x[1:])]; ks = np.arange(2, len(x) + 1)
    ax.semilogy(ks, dk, "-", marker=mk, color=col, lw=2, ms=5, mec="white", mew=1.0, label=lab)
    ax.text(ks[-1] + 0.15, dk[-1], f"ratio {dk[-1] / dk[-2]:.4f}", color=INK, fontsize=7.5, va="center")
ax.set_xlim(1.7, 8.6); ax.set_xlabel(r"level $k$"); ax.set_ylabel(r"$|\kappa_k-\kappa_{k-1}|$")
ax.legend(loc="upper right", fontsize=7.5)
save(fig, "hierarchy")
print("ok")
