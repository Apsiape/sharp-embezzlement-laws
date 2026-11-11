#!/usr/bin/env python3
"""Evaluate the XXZ corner-transfer-matrix spectrum in the pure Néel phase."""

import argparse
import bisect
import csv
import itertools
import math
import sys


A = math.pi**2 / 12


def distinct_partition_counts(maximum):
    """Use Euler's pentagonal recurrence, then multiply by (q²;q²)∞."""
    offsets = []
    for k in itertools.count(1):
        lower = k * (3 * k - 1) // 2
        if lower > maximum:
            break
        sign = 1 if k % 2 else -1
        offsets.extend(((lower, sign), (lower + k, sign)))

    ordinary = [1]
    distinct = [1]
    for level in range(1, maximum + 1):
        count = 0
        for offset, sign in offsets:
            if offset > level:
                break
            count += sign * ordinary[level - offset]
        ordinary.append(count)

        for offset, sign in offsets:
            if 2 * offset > level:
                break
            count -= sign * ordinary[level - 2 * offset]
        distinct.append(count)
    return distinct


def verify_small_counts(counts):
    maximum = min(256, len(counts) - 1)
    coefficients = [1] + [0] * maximum
    for part in range(1, maximum + 1):
        for level in range(maximum, part - 1, -1):
            coefficients[level] += coefficients[level - part]
    if coefficients != counts[: maximum + 1]:
        raise ArithmeticError("Distinct partition counts failed the direct check")


def partition_statistics(t):
    """Evaluate L(t) and L(t) - t L'(t) from the modular identity."""
    dual_spacing = 2 * math.pi**2 / t
    dual_levels = [
        (j, math.exp(-dual_spacing * j))
        for j in range(1, math.ceil(50 / dual_spacing) + 1)
    ]
    dual_log_partition = math.fsum(math.log1p(q) for _, q in dual_levels)
    dual_derivative = -math.fsum(j * q / (1 + q) for j, q in dual_levels)
    log_partition = A / t - math.log(2) / 2 + t / 24 - dual_log_partition
    entropy = (
        2 * A / t
        - math.log(2) / 2
        - dual_log_partition
        - dual_spacing * dual_derivative
    )
    return log_partition, entropy


def head_terms(ranks, t, d, log_partition):
    """Integrate the squared difference on the exact integer rank intervals."""
    log_dimension = math.log(d)
    source_level = 0
    target_level = 0
    start = 0
    while source_level < len(ranks):
        source_end = ranks[source_level]
        target_end = d * ranks[target_level]
        end = min(source_end, target_end)
        log_mass = math.log(end - start) - log_partition - t * source_level
        log_amplitude_ratio = (t * (source_level - target_level) - log_dimension) / 2
        yield math.exp(log_mass) * math.expm1(log_amplitude_ratio) ** 2
        if source_end == end:
            source_level += 1
        if target_end == end:
            target_level += 1
        start = end


def probability_tail_bound(level, t):
    """Upper bound for the Boltzmann probability of a level >= level."""
    if level == 0:
        return 1.0
    saddle = math.sqrt(A / level)
    if saddle >= t:
        return 1.0
    dual_spacing = 2 * math.pi**2 / t
    dual_log_partition_bound = math.exp(-dual_spacing) / -math.expm1(-dual_spacing)
    exponent = -t * (math.sqrt(level) - math.sqrt(A) / t) ** 2
    return math.exp(exponent + dual_log_partition_bound)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--etas", nargs="+", type=float, default=[0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.0025]
    )
    parser.add_argument("--targets", nargs="+", type=int, default=[2, 3, 5])
    args = parser.parse_args()
    if any(not math.isfinite(eta) or eta <= 0 for eta in args.etas):
        parser.error("all eta values must be finite and positive")
    if any(d < 2 for d in args.targets):
        parser.error("all target dimensions must be at least two")

    cutoffs = {
        eta: math.ceil((math.sqrt(A) / (2 * eta) + math.sqrt(44 / (2 * eta))) ** 2)
        for eta in args.etas
    }
    counts = distinct_partition_counts(max(cutoffs.values()))
    verify_small_counts(counts)
    cumulative = list(itertools.accumulate(counts))

    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(
        [
            "eta",
            "delta",
            "target_d",
            "entropy",
            "largest_probability",
            "embezzlement_error",
            "omitted_error_upper_bound",
            "normalized_error_entropy",
            "quadratic_residual",
            "cutoff_level",
        ]
    )
    for eta in args.etas:
        t = 2 * eta
        cutoff = cutoffs[eta]
        ranks = cumulative[: cutoff + 1]
        log_partition, entropy = partition_statistics(t)
        for d in args.targets:
            error = math.fsum(head_terms(ranks, t, d, log_partition))
            target_tail_level = bisect.bisect_right(ranks, ranks[-1] // d)
            tail = probability_tail_bound(cutoff + 1, t) + probability_tail_bound(
                target_tail_level, t
            )
            h_squared = math.log(d) ** 2
            leading_error = 3 * eta * h_squared / math.pi**2
            writer.writerow(
                [
                    eta,
                    math.cosh(eta),
                    d,
                    entropy,
                    math.exp(-log_partition),
                    error,
                    tail,
                    4 * error * entropy / h_squared,
                    (error - leading_error) / eta**2,
                    cutoff,
                ]
            )


if __name__ == "__main__":
    main()
