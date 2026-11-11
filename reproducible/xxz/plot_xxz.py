"""Reproduce the interacting XXZ spectrum figure from integer-rank data."""
from pathlib import Path
import csv
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--data", type=Path, default=Path(__file__).with_name("xxz_numerics.csv"))
parser.add_argument("--output", type=Path, default=Path(__file__).with_name("fig_xxz.pdf"))
args = parser.parse_args()
with args.data.open() as stream:
    rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
plt.rcParams.update({"font.family": "serif", "font.size": 8.3, "mathtext.fontset": "cm",
                     "axes.linewidth": .65, "legend.frameon": False, "legend.fontsize": 7.5,
                     "xtick.direction": "in", "ytick.direction": "in", "pdf.fonttype": 42})
fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.65), constrained_layout=True)
for d, color, marker in [(2, "#285f8f", "o"), (3, "#b55f30", "s"), (5, "#497e60", "^")]:
    selected = [row for row in rows if int(row["target_d"]) == d]
    x = np.array([1/row["eta"] for row in selected])
    y = np.array([row["embezzlement_error"] for row in selected])
    axes[0].loglog(x, y, marker+"-", color=color, lw=1, ms=3.5, label=fr"$d={d}$")
    axes[0].loglog(x, 3*np.log(d)**2/(np.pi**2*x), "--", color=color, lw=.75, alpha=.65)
    axes[1].semilogx(x, [row["normalized_error_entropy"] for row in selected], marker+"-",
                     color=color, lw=1, ms=3.5, label=fr"$d={d}$")
axes[0].set_ylabel(r"Embezzlement error $\mathcal{E}_d$")
axes[1].set_ylabel(r"$4S\mathcal{E}_d/(\ln d)^2$")
axes[1].axhline(1, color=".45", linestyle="--", lw=.8)
axes[1].set_ylim(.99, 1.27)
for i, ax in enumerate(axes):
    ax.set_xlabel(r"Inverse spacing $\eta^{-1}$, $\Delta=\cosh\eta$")
    ax.text(.02, .98, "("+chr(97+i)+")", transform=ax.transAxes, va="top")
    ax.legend(loc="upper right")
args.output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.output, metadata={"Title": "Exact XXZ near-critical embezzlement law"})
fig.savefig(args.output.with_suffix(".png"), dpi=180)
print(args.output)
