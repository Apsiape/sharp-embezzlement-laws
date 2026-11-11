import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/6fbe6fdd-66e6-555c-8896-876210039d7a/scratchpad/paper")
import json, glob, math, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT = "/home/claude/paper/"
plt.rcParams.update({"font.family": "serif", "font.size": 9, "axes.linewidth": 0.6, "mathtext.fontset": "cm",
                     "legend.frameon": False, "legend.fontsize": 7.6, "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
def load(pattern):
    rows = []
    for f in sorted(glob.glob(pattern)): rows += json.load(open(f))
    return rows
ising = sorted(load("data_ising_3_17.json"), key=lambda r: r["size"])
xx = sorted(load("data_xx_sq.json"), key=lambda r: r["size"])
floor = sorted(load("data_gapped_floor_*.json"), key=lambda r: -r["delta"])
gN = {d: sorted(load(f"data_gappedN_{d}.json"), key=lambda r: r["size"]) for d in (0.1, 0.01)}
cl = sorted(json.load(open("data_cl_fine.json")), key=lambda r: r["b"])
tw = json.load(open("data_towers.json")); ext = json.load(open("data_towers_ext.json"))
tw = {k: sorted(tw[k] + ext[k], key=lambda r: r["b"]) for k in ("1", "2")}
flat = 2 - math.sqrt(2); L2 = math.log(2)
c_i, c_x, c_g, ink, ink2 = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e"
def interp(rows):
    b = np.array([r["b"] for r in rows]); E = np.array([r["E"] for r in rows])
    return lambda x: float(np.exp(np.interp(math.log(x), np.log(b), np.log(E))))
clE, tw1, tw2 = interp(cl), interp(tw["1"]), interp(tw["2"])

# ---------------- Figure 1(b): embezzlement error versus system size (drawn at its printed size)
with plt.rc_context({"font.size": 8, "legend.fontsize": 7}):
    fig, ax = plt.subplots(figsize=(2.98, 1.95))
    ax.axhline(flat, color=ink, lw=0.9, ls=":")
    ax.text(8, flat + 0.014, r"product states and Bell pairs: $2-\sqrt{2}$", fontsize=7, color=ink)
    for d, mk, lab, mfc in ((0.1, "s", r"gapped chain, $\xi\approx 20$", "white"), (0.01, "D", r"gapped chain, $\xi\approx 200$", c_g)):
        ax.semilogx([r["size"] for r in gN[d]], [r["E"] for r in gN[d]], mk + "--", ms=2.8, lw=0.8, color=c_g, mfc=mfc, mew=0.8, label=lab)
    ax.semilogx([r["size"] for r in ising], [r["E"] for r in ising], "o-", ms=3.2, lw=1.0, color=c_i, label="critical Ising chain")
    ax.semilogx([r["size"] for r in xx], [r["E"] for r in xx], "^-", ms=3.4, lw=1.0, color=c_x, label="critical XX chain")
    ax.set_xlabel("number of sites $N$"); ax.set_ylabel(r"embezzlement error $\mathcal{E}$")
    ax.set_ylim(0, 0.68); ax.set_xlim(6, 4e5); ax.set_yticks([0, 0.2, 0.4, 0.6])
    h, l = ax.get_legend_handles_labels(); order = [2, 3, 0, 1]
    ax.legend([h[i] for i in order], [l[i] for i in order], loc="upper right", bbox_to_anchor=(1.0, 0.86), handlelength=2.2, labelspacing=0.3)
    ax.text(0.035, 0.05, "(b)", transform=ax.transAxes, fontsize=9)
    fig.tight_layout(pad=0.3); fig.savefig(OUT + "fig_main.pdf"); plt.close(fig)

# ---------------- Figure 2: the finite-size law
fig, axs = plt.subplots(1, 2, figsize=(6.7, 2.85))
ax = axs[0]
bb = np.exp(np.linspace(math.log(0.06), math.log(1.35), 300))
ax.plot(bb, [tw1(x) for x in bb], color=ink, lw=1.1, label="equally spaced levels, once")
ax.plot(bb, [tw2(x) for x in bb], color=ink, lw=1.1, ls="--", label="equally spaced levels, twice")
ax.plot(bb, [clE(x) for x in bb], color=ink2, lw=0.9, ls=":", label="smooth conformal spectrum")
ax.plot([r["b"] for r in floor], [r["E"] for r in floor], "s", ms=3.4, color=c_g, mfc="white", mew=0.9, label=r"gapped chains, $\xi=4,\dots,10^4$")
ax.plot([r["b"] for r in ising], [r["E"] for r in ising], "o", ms=3.2, color=c_i, label="critical Ising chain")
ax.plot([r["b"] for r in xx], [r["E"] for r in xx], "^", ms=3.4, color=c_x, label="critical XX chain")
ax.set_xlim(0, 1.3); ax.set_ylim(0, 0.5); ax.set_xlabel(r"$b=-\ln p_1$"); ax.set_ylabel(r"embezzlement error $\mathcal{E}$")
ax.legend(loc="upper right", handlelength=2.3); ax.text(0.035, 0.06, "(a)", transform=ax.transAxes, fontsize=10)
ax = axs[1]
f8 = lambda rows: ([r["b"] for r in rows], [8 * r["b"] * r["E"] / L2 ** 2 for r in rows])
ax.axhline(1.0, color=ink2, lw=0.6)
ax.semilogx(*f8(tw["1"]), color=ink, lw=1.1); ax.semilogx(*f8(tw["2"]), color=ink, lw=1.1, ls="--")
ax.semilogx(*f8(cl), color=ink2, lw=0.9, ls=":")
ax.semilogx(*f8(floor), "s", ms=3.4, color=c_g, mfc="white", mew=0.9)
ax.semilogx(*f8(ising), "o", ms=3.2, color=c_i); ax.semilogx(*f8(xx), "^", ms=3.4, color=c_x)
ax.set_xlim(0.1, 64); ax.set_ylim(0, 1.9); ax.set_xlabel(r"$b=-\ln p_1$"); ax.set_ylabel(r"$8b\,\mathcal{E}/(\ln 2)^2$")
ax.text(9, 1.07, "asymptotic law", fontsize=7.6, color=ink2); ax.text(0.035, 0.06, "(b)", transform=ax.transAxes, fontsize=10)
fig.tight_layout(pad=0.4, w_pad=1.4); fig.savefig(OUT + "fig_law.pdf"); plt.close(fig)

# ---------------- tables
def num(n): return f"{n:,}".replace(",", "{,}")
with open(OUT + "tab_chains.tex", "w") as f:
    for name, rows, sizes, twf in (("Critical Ising chain", ising, (16, 256, 4096, 65536, 131072), tw1), ("Critical XX chain", xx, (16, 256, 4096, 65536, 262144), tw2)):
        first = True
        for r in rows:
            if r["size"] not in sizes: continue
            f.write(f"{name if first else ''} & $N={num(r['size'])}$ & {r['S']:.3f} & {r['b']:.3f} & {r['E']:.4f} & {r['E']/twf(r['b']):.3f} & {r['E']/clE(r['b']):.2f} \\\\\n"); first = False
        f.write("\\addlinespace\n")
    first = True
    for r in floor:
        if r["delta"] not in (0.1, 0.01, 0.001, 0.0002): continue
        f.write(f"{'Gapped chain' if first else ''} & $\\xi={num(round(r['xi']))}$ & {r['S']:.3f} & {r['b']:.3f} & {r['E']:.4f} & {r['E']/tw1(r['b']):.3f} & {r['E']/clE(r['b']):.2f} \\\\\n"); first = False
pairs = json.load(open("data_pairs.json"))
with open(OUT + "tab_pairs.tex", "w") as f:
    for lam, lab in ((2 / 3, "$(2/3,\\,1/3)$ & 2"), (0.75, "$(3/4,\\,1/4)$ & 3"), (0.8, "$(4/5,\\,1/5)$ & 4")):
        vals = {d["n"]: d["E"] for d in pairs if abs(d["lam"] - lam) < 1e-9}
        f.write(lab + " & " + " & ".join(f"{vals[n]:.4f}" for n in (1, 16, 64, 256, 1000)) + " \\\\\n")
print("ratios to equally spaced prediction: ising", [(r["size"], round(r["E"] / tw1(r["b"]), 4)) for r in ising])
print("xx", [(r["size"], round(r["E"] / tw2(r["b"]), 4)) for r in xx])
print("gapped", [(round(r["xi"]), round(r["E"] / tw1(r["b"]), 4)) for r in floor])
print("ratio to smooth: ising", [(r["size"], round(r["E"] / clE(r["b"]), 3)) for r in ising][-4:], "xx", [(r["size"], round(r["E"] / clE(r["b"]), 3)) for r in xx][-6:], "gapped", [(round(r["xi"]), round(r["E"] / clE(r["b"]), 3)) for r in floor])
print("8bE/ln2^2: CL at b=1,2,3,4,8,16,64:", [round(8 * b * clE(b) / L2 ** 2, 4) for b in (1, 2, 3, 4, 8, 16, 64)], " towers at b=3.2:", round(8 * 3.2 * tw1(3.2) / L2 ** 2, 3), round(8 * 3.2 * tw2(3.2) / L2 ** 2, 3))
print("E*S: ising", [(r["size"], round(r["E"] * r["S"], 4)) for r in ising][-3:], "xx", [(r["size"], round(r["E"] * r["S"], 4)) for r in xx][-3:], "CL(b=0.57):", clE(0.57189) * 2 * 0.57189)
print("S vs ln N slope: ising", (ising[-1]["S"] - ising[-5]["S"]) / math.log(ising[-1]["size"] / ising[-5]["size"]), "xx", (xx[-1]["S"] - xx[-5]["S"]) / math.log(xx[-1]["size"] / xx[-5]["size"]))
print("uncertainty products E(4VarX+1): ising", [(r["size"], round(r["E"] * (4 * r["varX"] + 1), 3)) for r in ising][-2:], "xx", [(r["size"], round(r["E"] * (4 * r["varX"] + 1), 3)) for r in xx][-2:])
print("loss with Phi9: flat", 2 * flat / 243, "ising131072", 2 * ising[-1]["E"] / 243, "xx262144", 2 * xx[-1]["E"] / 243)
