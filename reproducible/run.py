"""Regenerate numerical data from the model parameters in the manuscript."""

import argparse
import json
import math
from pathlib import Path
from time import perf_counter

import numpy as np

from chains import ModeSpectrum, chain_row, ising_modes, row_from_modes
from conformal import pair_row, smooth_row, tower_rows


ROOT = Path(__file__).resolve().parent


def write_json(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rows, indent=2, allow_nan=False) + "\n")


def progress(row, started):
    print(f"{row['model']:10} N={row['size']:7} E={row['E']:.9f} "
          f"S={row['S']:.9f} mass={row['mass']:.11f} seconds={perf_counter() - started:.2f}",
          flush=True)


def run_chains(args):
    rows = []
    for power in range(args.min_power, args.max_power + 1):
        size = 2 ** power
        started = perf_counter()
        if args.model == "critical":
            if power == args.min_power:
                # Reproduce the XX point at the same smallest size as the Ising sweep.
                xx = chain_row("xx-doubled", size, threshold=args.threshold, count=args.modes)
                rows.append(xx)
                progress(xx, started)
                started = perf_counter()
            modes = ising_modes(size, args.modes or 60)
            doubled = ModeSpectrum(np.tile(modes.occupation, 2), 2 * modes.full_count)
            ising = row_from_modes(modes, "ising", size, args.threshold)
            xx = row_from_modes(doubled, "xx-doubled", 2 * size, args.threshold)
            rows.extend([ising, xx])
            progress(ising, started)
            progress(xx, started)
        else:
            row = chain_row(args.model, size, args.delta, args.threshold, args.modes)
            rows.append(row)
            progress(row, started)
        write_json(args.output, rows)


def run_gapped_floor(args):
    rows = []
    for delta in args.deltas:
        xi = 1 / math.asinh(delta / 2)
        power = math.ceil(math.log2(max(256, 8 * xi)))
        previous = None
        while power <= args.max_power:
            started = perf_counter()
            row = chain_row("gapped", 2 ** power, delta, args.threshold, args.modes)
            progress(row, started)
            difference = None if previous is None else abs(previous - row["E"])
            if power == args.max_power or (difference is not None and difference < 1e-6):
                row.update(change_on_doubling=difference,
                           within_doubling_tolerance=difference is not None and difference < 1e-6)
                rows.append(row)
                write_json(args.output, rows)
                break
            previous = row["E"]
            power += 1
        else:
            raise ValueError(f"No size allowed for delta={delta}. Increase --max-power.")


def run_theory(args):
    if not 0 < args.min_epsilon <= 0.25:
        raise ValueError("Use 0 < --min-epsilon <= 0.25.")
    legacy_fine = np.arange(0.02, 1.4001, 0.02).round(4).tolist()
    points = [1.5, 1.6, 1.7, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.2, 3.5, 4, 5,
              6, 8, 12, 16, 24, 32, 48, 64, 128]
    smooth = [smooth_row(b, args.dilation) for b in legacy_fine + points]
    towers = {}
    settings = ((1, 0.25, 70, (0.23, 0.21, 0.19, 0.175, 0.16, 0.148, 0.137, 0.127)),
                (2, 0.40, 60, (0.37, 0.34, 0.31, 0.285, 0.265, 0.25)))
    for copies, limit, grid_count, legacy_extension in settings:
        eps = np.geomspace(5., limit, grid_count).tolist()
        eps.extend(np.geomspace(limit, args.min_epsilon, 30).tolist()[1:])
        eps.extend(value for value in legacy_extension if value >= args.min_epsilon)
        eps = sorted(set(eps), reverse=True)
        towers[str(copies)] = list(tower_rows(eps, copies, args.dilation))
    write_json(args.output, {"smooth": smooth, "towers": towers})


def run_pairs(args):
    rows = [pair_row(p, n, args.dilation) for p in (2 / 3, 0.8, 0.9, 0.75)
            for n in args.counts]
    write_json(args.output, rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(required=True)
    chains = commands.add_parser("chains", help="Generate finite-chain data.")
    chains.add_argument("model", choices=("critical", "ising", "xx", "xx-doubled", "gapped"))
    chains.add_argument("--min-power", type=int, default=3)
    chains.add_argument("--max-power", type=int, default=10)
    chains.add_argument("--delta", type=float, default=0.1)
    chains.set_defaults(function=run_chains)
    floor = commands.add_parser("gapped-floor", help="Estimate the large-size gapped values.")
    floor.add_argument("deltas", type=float, nargs="+")
    floor.add_argument("--max-power", type=int, default=18)
    floor.set_defaults(function=run_gapped_floor)
    for command in (chains, floor):
        command.add_argument("--modes", type=int, help="Number of one-body Lanczos modes.")
        command.add_argument("--threshold", type=float, default=1e-11)
        command.add_argument("--output", type=Path, required=True)
    theory = commands.add_parser("theory", help="Regenerate the conformal models.")
    theory.add_argument("--min-epsilon", type=float, default=0.02)
    theory.set_defaults(function=run_theory)
    pairs = commands.add_parser("pairs", help="Regenerate repeated-pair data.")
    pairs.add_argument("--counts", type=int, nargs="+", default=(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1000))
    pairs.set_defaults(function=run_pairs)
    for command in (theory, pairs):
        command.add_argument("--dilation", type=int, default=2)
        command.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
