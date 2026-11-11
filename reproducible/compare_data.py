"""Report agreement with the archived numerical supplement."""

import json
import platform
from pathlib import Path

import matplotlib
import numpy
import scipy


ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT / "original" / "supplement" / "data"
GENERATED = ROOT / "generated"


def load(path):
    return json.loads(path.read_text())


def compare(rows, reference, key, fields):
    target = {key(row): row for row in rows}
    matches = [(target[key(row)], row) for row in reference]
    return dict(rows=len(matches), max_abs_difference={
        field: max(abs(a[field] - b[field]) for a, b in matches) for field in fields
    })


def comparison():
    critical = load(GENERATED / "critical.json")
    theory = load(GENERATED / "theory.json")
    chain_fields = ("S", "b", "E", "varX", "varK", "mass")
    report = {"runtime": {"python": platform.python_version(), "numpy": numpy.__version__,
                          "scipy": scipy.__version__, "matplotlib": matplotlib.__version__}}
    for model, filename in (("ising", "data_ising_3_17.json"), ("xx-doubled", "data_xx_sq.json")):
        report[model] = compare([row for row in critical if row["model"] == model],
                                load(ARCHIVE / filename), lambda row: row["size"], chain_fields)
    report["xx_against_archived_direct"] = compare([row for row in critical if row["model"] == "xx-doubled"],
                                                    load(ARCHIVE / "data_xx_3_16.json"),
                                                    lambda row: row["size"], chain_fields)

    floor = load(GENERATED / "gapped-floor-samples.json")
    deltas = {row["delta"] for row in floor}
    originals = [row for path in ARCHIVE.glob("data_gapped_floor_*.json") for row in load(path)
                 if row["delta"] in deltas]
    report["gapped_floor_samples"] = compare(floor, originals, lambda row: row["delta"], chain_fields)
    report["pairs"] = compare(load(GENERATED / "pairs.json"), load(ARCHIVE / "data_pairs.json"),
                              lambda row: (round(row["lam"], 12), row["n"]), ("E", "S"))
    report["smooth"] = compare(theory["smooth"], load(ARCHIVE / "data_cl_fine.json"),
                               lambda row: round(row["b"], 12), ("E", "S", "mass"))
    for copies in ("1", "2"):
        archived = load(ARCHIVE / "data_towers.json")[copies] + load(ARCHIVE / "data_towers_ext.json")[copies]
        report[f"tower_{copies}"] = compare(theory["towers"][copies], archived,
                                           lambda row: round(row["eps"], 12), ("E", "S", "b", "mass"))
    return report


if __name__ == "__main__":
    print(json.dumps(comparison(), indent=2))
