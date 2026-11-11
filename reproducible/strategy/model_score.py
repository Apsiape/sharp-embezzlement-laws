"""Certified uniform-game score for a doubly repeated geometric mode list."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from functools import cache
import json
from math import isqrt
from pathlib import Path


PLACES = 45
SCALE = 10 ** PLACES


def ceil_div(numerator: int, denominator: int) -> int:
    return -(-numerator // denominator)


@dataclass(frozen=True, slots=True)
class Interval:
    """Endpoints stored as integer multiples of 10**(-PLACES)."""

    low: int
    high: int

    @classmethod
    def ratio(cls, numerator: int, denominator: int = 1) -> Interval:
        return cls(numerator * SCALE // denominator,
                   ceil_div(numerator * SCALE, denominator))

    def __add__(self, other: Interval) -> Interval:
        return Interval(self.low + other.low, self.high + other.high)

    def __sub__(self, other: Interval) -> Interval:
        return Interval(self.low - other.high, self.high - other.low)

    def __mul__(self, other: Interval) -> Interval:
        """Multiply nonnegative intervals."""
        assert self.low >= 0 and other.low >= 0
        return Interval(self.low * other.low // SCALE,
                        ceil_div(self.high * other.high, SCALE))

    def times(self, integer: int) -> Interval:
        return Interval(self.low * integer, self.high * integer)

    def divide(self, integer: int) -> Interval:
        return Interval(self.low // integer, ceil_div(self.high, integer))

    def reciprocal(self) -> Interval:
        assert self.low > 0
        return Interval(SCALE * SCALE // self.high,
                        ceil_div(SCALE * SCALE, self.low))

    def square_root(self) -> Interval:
        assert self.low >= 0
        lower = isqrt(self.low * SCALE)
        upper = isqrt(self.high * SCALE)
        return Interval(lower, upper + (upper * upper < self.high * SCALE))

    def nonnegative(self) -> Interval:
        assert self.high >= 0
        return Interval(max(0, self.low), self.high)

    def report(self) -> dict[str, str]:
        def display(value: int) -> str:
            whole, fractional = divmod(value, SCALE)
            return f"{whole}.{fractional:0{PLACES}d}"

        assert self.low >= 0
        return {"lower": display(self.low), "upper": display(self.high)}


ZERO = Interval.ratio(0)
ONE = Interval.ratio(1)
SQRT_TWO = Interval.ratio(2).square_root()
BASE = Interval.ratio(32767, 32768)


@cache
def amplitude(exponent: int) -> Interval:
    """Enclose (32767/32768)**exponent by directed integer arithmetic."""
    if exponent == 0:
        return ONE
    half = amplitude(exponent // 2)
    square = half * half
    return square * BASE if exponent % 2 else square


def rank_counts(exponents: list[int], cutoff: int) -> list[int]:
    """Coefficients of the truncated product of (1 + x**k)**2."""
    counts = [1] + [0] * cutoff
    top = 0
    for exponent in sorted(exponents):
        if exponent == 0:
            counts[0] *= 4
            continue
        if exponent > cutoff:
            break
        doubled = 2 * exponent
        top = min(cutoff, top + doubled)
        for energy in range(top, doubled - 1, -1):
            counts[energy] += 2 * counts[energy - exponent] + counts[energy - doubled]
        for energy in range(min(top, doubled - 1), exponent - 1, -1):
            counts[energy] += 2 * counts[energy - exponent]
    return counts


def run_statistics(counts: list[int]) -> tuple[int, Interval, Interval, Interval]:
    """Return rank, mass, duplication overlap and nine-block variance before normalization."""
    runs = []
    rank = previous = pending = 0
    coefficient = ONE
    mass = variance = total = squares = ZERO
    for energy, count in enumerate(counts):
        if not count:
            continue
        coefficient = coefficient * amplitude(energy - previous)
        previous = energy
        probability = coefficient * coefficient
        rank += count
        runs.append((rank, coefficient))
        mass += probability.times(count)

        # Complete the pending mixed block. Full constant blocks have variance zero.
        if pending:
            length = min(9 - pending, count)
            total += coefficient.times(length)
            squares += probability.times(length)
            pending += length
            count -= length
            if pending == 9:
                variance += (squares - (total * total).divide(9)).nonnegative()
                pending = 0
        if count:
            assert pending == 0
            pending = count % 9
            total = coefficient.times(pending)
            squares = probability.times(pending)
    if pending:
        variance += (squares - (total * total).divide(9)).nonnegative()

    overlap = ZERO
    start = direct = dilated = 0
    while direct < len(runs):
        direct_end, direct_value = runs[direct]
        dilated_end, dilated_value = runs[dilated]
        dilated_end *= 2
        end = min(direct_end, dilated_end)
        overlap += (direct_value * dilated_value).times(end - start)
        direct += end == direct_end
        dilated += end == dilated_end
        start = end
    return rank, mass, overlap, variance


def score(exponents: list[int], cutoff: int, distance_bound: Fraction = Fraction(0)) -> dict:
    """Enclose the prefix statistics and the strategy's loss upper bound."""
    if cutoff < 0 or any(type(k) is not int or k < 0 for k in exponents):
        raise ValueError("Use nonnegative integer exponents and cutoff.")
    if distance_bound < 0:
        raise ValueError("Use a nonnegative vector-distance bound.")
    inverse_peak = ONE
    for exponent in exponents:
        term = ONE + amplitude(2 * exponent)
        inverse_peak = inverse_peak * term * term
    peak = inverse_peak.reciprocal()

    rank, mass, overlap, variance = run_statistics(rank_counts(exponents, cutoff))
    # The dilation has the same mass and its overlap has the factor 1/sqrt(2).
    error = peak * (mass.times(2) - SQRT_TWO * overlap).nonnegative()
    variance = peak * variance
    tail = (ONE - peak * mass).nonnegative()

    # Homogeneity makes normalization of the prefix unnecessary.
    combined = error.square_root() + variance.times(8).divide(9).square_root()
    raw_upper = (combined * combined).times(2) + variance.times(72) + tail.times(9)
    loss_upper = raw_upper.divide(243)
    raw_floor = Interval.ratio(3) - SQRT_TWO.times(2)
    floor = raw_floor.divide(243)
    distance = Interval.ratio(distance_bound.numerator, distance_bound.denominator)
    transferred_amplitude = loss_upper.square_root() + distance
    transferred_loss = transferred_amplitude * transferred_amplitude
    return {
        "exponents": exponents,
        "cutoff": cutoff,
        "rank": rank,
        "peak": peak.report(),
        "tail": tail.report(),
        "prefix_duplication_error": error.report(),
        "prefix_nine_block_variance": variance.report(),
        "raw_loss_bound": raw_upper.report(),
        "loss_bound": loss_upper.report(),
        "uniform_tracial_floor": floor.report(),
        "model_separates": raw_upper.high < raw_floor.low,
        "distance_bound": str(distance_bound),
        "transferred_loss_bound": transferred_loss.report(),
        "distance_bound_certifies_separation": transferred_loss.high < floor.low,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSON exponent list or object with exponents")
    parser.add_argument("--cutoff", type=int, default=400_000)
    parser.add_argument("--distance-bound", type=Fraction, default=Fraction(0),
                        help="certified upper bound on aligned vector distance")
    parser.add_argument("--require-separation", action="store_true",
                        help="fail unless the model and transferred bounds both beat the flat-state floor")
    args = parser.parse_args()
    if not __debug__:
        raise SystemExit("Run the certificate without Python optimization (-O or PYTHONOPTIMIZE).")
    values = json.loads(args.input.read_text())
    exponents = values["exponents"] if isinstance(values, dict) else values
    result = score(exponents, args.cutoff, args.distance_bound)
    print(json.dumps(result, indent=2))
    if args.require_separation and not (result["model_separates"] and
                                        result["distance_bound_certifies_separation"]):
        raise SystemExit("The supplied certificate does not prove separation.")


if __name__ == "__main__":
    main()
