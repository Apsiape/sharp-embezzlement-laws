"""Regenerate the manuscript figures, tables, and plotting data."""

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from game import CLASSICAL_LOSS
from spectrum import ASSISTED_COEFFICIENT, FLAT_ERROR, LOG2


ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT / "original" / "supplement" / "data"
GENERATED = ROOT / "generated"
ISING, XX, GAPPED, GREY = "#1c67a7", "#b44b2d", "#25846a", "#494e53"
LABELS = {"ising": "Critical Ising", "xx": "Critical XX", "gapped": "Staggered XX"}


def read(path):
    return json.loads(path.read_text())


def ordered(rows, key):
    return sorted(rows, key=lambda row: row[key])


def points(rows, key):
    return [row[key] for row in rows]


def load_data(archive_only):
    if archive_only:
        critical = {"ising": read(ARCHIVE / "data_ising_3_17.json"),
                    "xx": read(ARCHIVE / "data_xx_sq.json")}
        towers = read(ARCHIVE / "data_towers.json")
        for copies, extension in read(ARCHIVE / "data_towers_ext.json").items():
            towers[copies] = ordered(towers[copies] + extension, "b")
        smooth = read(ARCHIVE / "data_cl_fine.json")
        pairs = read(ARCHIVE / "data_pairs.json")
    else:
        rows = read(GENERATED / "critical.json")
        critical = {"ising": [row for row in rows if row["model"] == "ising"],
                    "xx": [row for row in rows if row["model"] == "xx-doubled"]}
        theory = read(GENERATED / "theory.json")
        towers, smooth = theory["towers"], theory["smooth"]
        pairs = read(GENERATED / "pairs.json")
    floor = [row for path in ARCHIVE.glob("data_gapped_floor_*.json") for row in read(path)]
    gapped = {delta: read(ARCHIVE / f"data_gappedN_{delta}.json") for delta in (0.1, 0.01)}
    if not archive_only:
        regenerated_floor = GENERATED / "gapped-floor.json"
        if regenerated_floor.is_file():
            floor = read(regenerated_floor)
        for delta in gapped:
            regenerated = GENERATED / f"gapped-{delta}.json"
            if regenerated.is_file():
                gapped[delta] = read(regenerated)
    return critical, ordered(floor, "b"), gapped, smooth, towers, pairs


def save(fig, output):
    fig.savefig(output, bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)


def make_size_figure(output, critical, gapped):
    fig, ax = plt.subplots(figsize=(3.22, 2.18), layout="constrained")
    ax.axhline(FLAT_ERROR, color=GREY, lw=0.85, ls=":")
    ax.text(9, FLAT_ERROR + 0.019, r"Flat spectrum: $2-\sqrt{2}$", color=GREY)
    for model, color, marker in (("ising", ISING, "o"), ("xx", XX, "^")):
        rows = critical[model]
        ax.semilogx(points(rows, "size"), points(rows, "E"), marker + "-",
                    color=color, label=LABELS[model], ms=3.1, lw=1.1)
    for delta, marker, label in ((0.1, "s", r"Gapped, $\xi\approx20$"),
                                 (0.01, "D", r"Gapped, $\xi\approx200$")):
        rows = gapped[delta]
        ax.semilogx(points(rows, "size"), points(rows, "E"), marker + "--", color=GAPPED,
                    label=label, ms=2.8, lw=0.9, mfc="white" if delta == 0.1 else GAPPED)
    ax.set(xlabel=r"Chain length $N$", ylabel=r"Embezzlement error $\mathcal{E}$",
           xlim=(6, 4e5), ylim=(0, 0.68), yticks=(0, 0.2, 0.4, 0.6))
    ax.legend(loc="upper right", bbox_to_anchor=(1, 0.84), ncols=1)
    ax.text(0.04, 0.05, "(b)", transform=ax.transAxes)
    save(fig, output)


def make_law_figure(output, critical, floor, smooth, towers):
    fig, axes = plt.subplots(1, 2, figsize=(6.75, 2.55), layout="constrained")
    theory = ((towers["1"], "-", "Single tower"), (towers["2"], "--", "Double tower"),
              (smooth, ":", "Smooth spectrum"))
    samples = ((critical["ising"], ISING, "o", LABELS["ising"]),
               (critical["xx"], XX, "^", LABELS["xx"]), (floor, GAPPED, "s", "Gapped limit estimates"))
    for index, ax in enumerate(axes):
        for rows, style, label in theory:
            values = [row["E"] * (8 * row["b"] / LOG2 ** 2 if index else 1) for row in rows]
            ax.plot(points(rows, "b"), values, style, color=GREY, label=label, lw=1.15)
        for rows, color, marker, label in samples:
            values = [row["E"] * (8 * row["b"] / LOG2 ** 2 if index else 1) for row in rows]
            ax.plot(points(rows, "b"), values, marker, color=color, label=label, ms=3.2)
        ax.set_xlabel(r"Largest probability parameter $b=-\ln p_1$")
        ax.text(0.035, 0.06, ("(a)", "(b)")[index], transform=ax.transAxes)
    axes[0].set(xlim=(0, 1.3), ylim=(0, 0.53), ylabel=r"Embezzlement error $\mathcal{E}$")
    axes[0].legend(loc="upper right", ncols=1, handlelength=2.)
    axes[1].set(xscale="log", xlim=(0.1, 64), ylim=(0, 1.9), ylabel=r"$8b\,\mathcal{E}/(\ln2)^2$")
    axes[1].axhline(1, color="#898d91", lw=0.7, zorder=0)
    axes[1].text(3.6, 0.86, "Leading asymptotic", color=GREY)
    save(fig, output)


def make_diagnostics(output, critical, smooth, towers):
    fig, (ax, loss) = plt.subplots(1, 2, figsize=(6.75, 2.65), layout="constrained")
    for rows, style, label in ((towers["1"], "-", "Single tower"),
                               (towers["2"], "--", "Double tower"),
                               (smooth, ":", "Smooth spectrum")):
        ax.semilogx(points(rows, "S"), [4 * row["S"] * row["E"] / LOG2 ** 2 for row in rows],
                    style, label=label, color=GREY, lw=1.1)
    ax.axhline(1, color="#8a9094", lw=0.65, zorder=0)
    for model, color, marker in (("ising", ISING, "o"), ("xx", XX, "^")):
        rows = critical[model]
        ax.semilogx(points(rows, "S"), [4 * row["S"] * row["E"] / LOG2 ** 2 for row in rows],
                    marker, label=LABELS[model], color=color, ms=3.1)
        loss.loglog(points(rows, "size"), points(rows, "loss_explicit"), marker + "-",
                    color=color, ms=2.5, lw=1.1)
        loss.loglog(points(rows, "size"), [ASSISTED_COEFFICIENT * row["E"] for row in rows],
                    marker + "--", color=color, ms=2.5, lw=1.)
    ax.set(xlabel=r"Entanglement entropy $S$ (nats)", ylabel=r"$4S\,\mathcal{E}/(\ln2)^2$",
           xlim=(0.26, 160), ylim=(0.7, 2.12))
    ax.legend(loc="upper right", handlelength=2., fontsize=7)
    ax.text(0.94, 0.07, "(a)", transform=ax.transAxes)
    loss.axhline(CLASSICAL_LOSS, color="#8a9094", lw=0.8, ls=":", zorder=0)
    loss.set(xlabel=r"Chain length $N$", ylabel="Strategy loss", xlim=(6, 4e5), ylim=(5e-4, 0.29))
    handles = [Line2D([0], [0], color=GREY, lw=1.1, label=r"Explicit, no shared $\Phi_9$"),
               Line2D([0], [0], color=GREY, lw=1.1, ls="--", label=r"With shared $\Phi_9$"),
               Line2D([0], [0], color="#8a9094", lw=0.8, ls=":", label="Classical optimum")]
    loss.legend(handles=handles, loc="center left", handlelength=2., fontsize=7)
    loss.text(0.04, 0.07, "(b)", transform=loss.transAxes)
    save(fig, output)


def export_data(output, critical, floor, gapped, smooth, towers):
    groups = list(critical.items()) + [("gapped floor", floor), ("smooth", smooth)]
    groups += [(f"gapped delta={delta}", rows) for delta, rows in gapped.items()]
    groups += [(f"tower {copies}", rows) for copies, rows in towers.items()]
    columns = ["series", "size", "delta", "xi", "eps", "S", "b", "E", "ratio_8bE", "ratio_4SE",
               "loss_explicit", "loss_assisted", "loss_B9", "loss_upper_prefix", "E_lower_with_modes",
               "E_upper_with_modes", "E_lower", "E_upper", "mass", "tail_mass_with_modes_bound"]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for name, rows in groups:
            for row in rows:
                entry = dict(row, series=name, ratio_8bE=8 * row["b"] * row["E"] / LOG2 ** 2,
                             ratio_4SE=4 * row["S"] * row["E"] / LOG2 ** 2)
                if "loss_explicit" in entry:
                    entry["loss_assisted"] = ASSISTED_COEFFICIENT * row["E"]
                writer.writerow(entry)


def interpolate(b, rows):
    return math.exp(float(np.interp(math.log(b), np.log(points(rows, "b")), np.log(points(rows, "E")))))


def export_tables(output, critical, floor, smooth, towers, pairs):
    chain_lines = []
    selections = {"ising": (16, 256, 4096, 65536, 131072), "xx": (16, 256, 4096, 65536, 262144)}
    for model, copies in (("ising", "1"), ("xx", "2"), ("gapped", "1")):
        selected = ([row for row in critical[model] if row["size"] in selections[model]]
                    if model != "gapped" else [row for row in floor if row["delta"] in (0.1, 0.01, 0.001, 0.0002)])
        for i, row in enumerate(selected):
            name = LABELS[model] if i == 0 else ""
            size = f"{row['size']:,}" if model != "gapped" else f"{round(row['xi']):,}"
            label = ("N=" if model != "gapped" else r"\xi=") + size.replace(",", "{,}")
            pred = interpolate(row["b"], towers[copies])
            cl = interpolate(row["b"], smooth)
            chain_lines.append(f"{name} & ${label}$ & {row['S']:.3f} & {row['b']:.3f} & "
                               f"{row['E']:.4f} & {row['E']/pred:.3f} & {row['E']/cl:.2f} " + r"\\")
        chain_lines.append(r"\addlinespace")
    (output / "tab_chains.tex").write_text("\n".join(chain_lines[:-1]) + "\n")

    pair_lines = []
    for weight, ratio in ((2 / 3, 2), (3 / 4, 3), (4 / 5, 4)):
        by_count = {row["n"]: row["E"] for row in pairs if abs(row["lam"] - weight) < 1e-10}
        numbers = " & ".join(f"{by_count[n]:.4f}" for n in (1, 16, 64, 256, 1000))
        pair_lines.append(f"$({ratio}/{ratio + 1},1/{ratio + 1})$ & {ratio} & {numbers} " + r"\\")
    (output / "tab_pairs.tex").write_text("\n".join(pair_lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=GENERATED / "figures")
    parser.add_argument("--archive-only", action="store_true", help="Use the supplied legacy data only.")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "serif", "font.size": 8.3, "mathtext.fontset": "cm",
                         "axes.linewidth": 0.65, "legend.frameon": False, "legend.fontsize": 7.1,
                         "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
                         "ytick.right": True, "pdf.fonttype": 42})
    critical, floor, gapped, smooth, towers, pairs = load_data(args.archive_only)
    make_size_figure(args.output / "fig_main.pdf", critical, gapped)
    make_law_figure(args.output / "fig_law.pdf", critical, floor, smooth, towers)
    if not args.archive_only:
        make_diagnostics(args.output / "fig_diagnostics.pdf", critical, smooth, towers)
    export_data(args.output / "figure_data.csv", critical, floor, gapped, smooth, towers)
    export_tables(args.output, critical, floor, smooth, towers, pairs)


if __name__ == "__main__":
    main()
