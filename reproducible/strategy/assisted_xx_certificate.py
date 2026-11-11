"""Certify the assisted uniform-game loss of open XX chains at powers 18–20.

The shared state is the middle-cut XX ground state tensored with Phi_9.
Floating-point eigensolvers only propose roots and eigenpairs for the
certificate. The separate nominal values are floating-point estimates.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh, eigh_tridiagonal

from certify_xx84_modes import arctan_enclosure
from model_score import BASE, Interval, ONE, SCALE, SQRT_TWO, amplitude, ceil_div, rank_counts, run_statistics

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from spectrum import error as duplication_error, mode_probabilities


ORDER = 24
ROOT_SCALE = 10 ** 28
ROOT_SEARCH_RADIUS = ROOT_SCALE // 10 ** 12
ENERGY_CUTOFF = 400_000
NOMINAL_RELATIVE_THRESHOLD = 1e-14
EXACT_FLOOR = (Interval.ratio(3) - SQRT_TWO.times(2)).divide(243)
REPORTED_LOSS_BOUNDS = {18: 746, 19: 703, 20: 663}
SCORE_DISTANCE_BUDGETS = {
    18: Interval.ratio(1363, 10 ** 8).high,
    19: Interval.ratio(9187, 10 ** 9).high,
    20: Interval.ratio(1684, 10 ** 8).high,
}
ERROR_CAPS = {
    "covariance_compression_singular_l2_upper": Interval.ratio(795, 10 ** 10).low,
    "kernel_enclosure_frobenius_upper": Interval.ratio(332, 10 ** 20).low,
    "omitted_squared_singular_mass_upper": Interval.ratio(3, 10 ** 15).low,
    "retained_squared_singular_error_upper": Interval.ratio(121, 10 ** 13).low,
    "total_one_parity_singular_l2_upper": Interval.ratio(355, 10 ** 8).low,
    "omitted_geometric_probability_upper": Interval.ratio(130, 10 ** 10).low,
}

# Fixed model proposals. Every mode occurs twice, with amplitude ratio
# (32767/32768)**k. The spectral residual check verifies the proposals.
EXPONENTS = {
    18: [11037, 33596, 57136, 81569, 106688, 132338, 158423, 184879,
         211661, 238735, 266074, 293658, 321469, 349490, 377711, 406118,
         434702, 463456, 492369, 521437, 550651, 580008, 609500, 639124],
    19: [10537, 32018, 54344, 77460, 101191, 125398, 149994, 174920,
         200135, 225608, 251316, 277239, 303362, 329670, 356152, 382798,
         409599, 436548, 463637, 490860, 518212, 545687, 573281, 600990],
    20: [10080, 30585, 51822, 73761, 96253, 119175, 142448, 166016,
         189843, 213900, 238165, 262622, 287256, 312054, 337007, 362104,
         387338, 412703, 438192, 463798, 489518, 515347, 541280, 567314],
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def signed_product(left: Interval, right: Interval) -> Interval:
    products = (left.low * right.low, left.low * right.high,
                left.high * right.low, left.high * right.high)
    return Interval(min(products) // SCALE, ceil_div(max(products), SCALE))


def absolute(interval: Interval) -> Interval:
    return Interval(max(0, interval.low, -interval.high),
                    max(-interval.low, interval.high))


def ceil_sqrt(value: int) -> int:
    root = math.isqrt(value)
    return root + (root * root < value)


def as_interval(value: Fraction) -> Interval:
    return Interval.ratio(value.numerator, value.denominator)


def upper_text(value: int) -> str:
    return Interval(0, value).report()["upper"]


def determinant_coefficients(squares: list[Fraction]) -> list[int]:
    """Integer coefficients for Horner evaluation of p_k(n / ROOT_SCALE)."""
    previous = [Fraction(1)]
    current = [Fraction(0), Fraction(1)]
    for square in squares:
        following = [Fraction(0)] + current
        for index, coefficient in enumerate(previous):
            following[index] -= square * coefficient
        previous, current = current, following
    denominator = math.lcm(*(coefficient.denominator for coefficient in current))
    return [int(coefficient * denominator) * ROOT_SCALE ** index
            for index, coefficient in enumerate(reversed(current))]


def determinant_at(point: int, coefficients: list[int]) -> int:
    value = 0
    for coefficient in coefficients:
        value = value * point + coefficient
    return value


def root_interval(proposal: float, coefficients: list[int]) -> tuple[int, int]:
    center = round(proposal * ROOT_SCALE)
    low, high = center - ROOT_SEARCH_RADIUS, center + ROOT_SEARCH_RADIUS
    low_value = determinant_at(low, coefficients)
    high_value = determinant_at(high, coefficients)
    require(low_value * high_value < 0, "Jacobi root proposal has no sign-changing bracket")
    positive_at_low = low_value > 0
    while high - low > 1:
        middle = (low + high) // 2
        value = determinant_at(middle, coefficients)
        if value == 0:
            return middle, middle
        if (value > 0) == positive_at_low:
            low = middle
        else:
            high = middle
    return low, high


def quadrature_bin(first: int, length: int) -> tuple[list[Interval], list[Interval]]:
    """Enclose the nodes and weights for counting measure on a single bin."""
    if length <= ORDER:
        return [Interval.ratio(i) for i in range(first, first + length)], [ONE] * length

    squares = [Fraction(j ** 2 * (length ** 2 - j ** 2),
                        4 * (4 * j ** 2 - 1) * length ** 2) for j in range(1, ORDER)]
    guesses = eigh_tridiagonal(np.zeros(ORDER), np.sqrt([float(x) for x in squares]),
                               eigvals_only=True)
    coefficients = determinant_coefficients(squares)
    brackets = [root_interval(float(guess), coefficients) for guess in guesses]
    require(len(brackets) == ORDER and all(left[1] < right[0]
            for left, right in zip(brackets, brackets[1:])), "Jacobi root brackets overlap")

    nodes, weights = [], []
    beta_squares = [as_interval(square) for square in squares]
    center = Interval.ratio(2 * first + length - 1, 2)
    for low, high in brackets:
        step = SCALE // ROOT_SCALE
        root = Interval(low * step, high * step)
        node = center + root.times(length)
        require(first * SCALE <= node.low <= node.high <= (first + length - 1) * SCALE,
                "Jacobi node enclosure leaves its bin")
        nodes.append(node)

        # With first component one, the j-th component is the leading
        # determinant p_j(root) divided by beta_1 * ... * beta_j.
        previous, current = ONE, root
        denominator = norm_squared = ONE
        for beta_square in beta_squares:
            denominator *= beta_square
            magnitude = absolute(current)
            norm_squared += magnitude * magnitude * denominator.reciprocal()
            previous, current = (current, signed_product(root, current)
                                 - signed_product(beta_square, previous))
        weights.append(norm_squared.reciprocal().times(length))
    return nodes, weights


def pi_interval() -> Interval:
    low5, high5 = arctan_enclosure(5)
    low239, high239 = arctan_enclosure(239)
    return Interval(as_interval(16 * low5 - 4 * high239).low,
                    as_interval(16 * high5 - 4 * low239).high)


TAYLOR_REMAINDER = Interval.ratio(1, 10 ** 18).high
TAYLOR_COEFFICIENTS = {
    shift: [Interval.ratio((-1) ** k, math.factorial(2 * k + shift))
            for k in range(12, -1, -1)]
    for shift in (0, 1)
}


def taylor_even(argument: Interval, shift: int) -> Interval:
    """Enclose cosine (shift 0) or sinc (shift 1) on arguments in [0, 2]."""
    require(0 <= argument.low <= argument.high <= 2 * SCALE, "Taylor argument exceeds its range")
    square = argument * argument
    coefficients = TAYLOR_COEFFICIENTS[shift]
    value = coefficients[0]
    for coefficient in coefficients[1:]:
        value = signed_product(value, square) + coefficient
    return Interval(value.low - TAYLOR_REMAINDER, value.high + TAYLOR_REMAINDER)


def enclosed_kernel(size: int) -> tuple[list[list[int]], int]:
    """Return the fixed-point center matrix and a Frobenius error upper bound."""
    require(SCALE % ROOT_SCALE == 0, "the interval grid cannot represent the root brackets")
    require(2 ** 26 * 10 ** 18 < math.factorial(26), "Taylor remainder exceeds its allowance")
    nodes, weights = [], []
    first = 1
    while first <= size // 4:
        length = min(first, size // 4 - first + 1)
        bin_nodes, bin_weights = quadrature_bin(first, length)
        nodes.extend(bin_nodes)
        weights.extend(bin_weights)
        first += length

    dimension = len(nodes)
    centers = [[0] * dimension for _ in nodes]
    squared_radius = 0
    denominator = 2 * (size + 1)
    frequency = pi_interval().divide(denominator)
    twice_nodes = [node.times(2) for node in nodes]
    weight_roots = [weight.square_root() for weight in weights]
    for i, (left, left_weight) in enumerate(zip(twice_nodes, weight_roots, strict=True)):
        for j in range(i + 1):
            right = twice_nodes[j]
            sine_argument = frequency * (left + right - ONE)
            cosine_argument = frequency * absolute(left - right)
            sine = sine_argument * taylor_even(sine_argument, 1)
            cosine = taylor_even(cosine_argument, 0)
            kernel = (sine.reciprocal() + cosine.reciprocal()).divide(denominator)
            entry = kernel * left_weight * weight_roots[j]
            center = (entry.low + entry.high) // 2
            centers[i][j] = centers[j][i] = center
            squared_radius += (1 if i == j else 2) * (entry.high - center) ** 2
    return centers, ceil_sqrt(squared_radius)


def compression_bound(size: int) -> Interval:
    """Analytic dyadic polynomial error for the exact weighted quadrature."""
    # The reciprocal kernel has Frobenius norm less than one on the
    # grid: (1 + log(N/2)) / (4*pi**2) < 20 / 36 < 1 for N <= 2**20.
    require(4 <= size <= 2 ** 20 and size % 4 == 0, "chain outside compression bound")
    require(pi_interval().low > 3 * SCALE, "pi enclosure cannot bound the reciprocal kernel")
    return Interval.ratio(15, 4 * 3 ** ORDER) + Interval.ratio(4, 3 * 2 ** ORDER)


def nearest_integer(value: float) -> int:
    """Round a binary64 proposal to the fixed-point grid."""
    numerator, denominator = float(value).as_integer_ratio()
    return (2 * numerator * SCALE + denominator) // (2 * denominator)


def checked_eigenvalue(matrix: list[list[int]], vector: np.ndarray) -> Interval:
    """Enclose an eigenvalue using the rounded proposal's exact Rayleigh quotient."""
    entries = [nearest_integer(entry) for entry in vector]
    norm_squared = sum(entry ** 2 for entry in entries)
    require(norm_squared > 0, "zero eigenvector proposal")
    applied = [sum(a * b for a, b in zip(row, entries, strict=True)) for row in matrix]
    rayleigh_numerator = sum(entry * result for entry, result in zip(entries, applied, strict=True))
    center = (2 * rayleigh_numerator + norm_squared) // (2 * norm_squared)
    residual_squared = sum((result - center * entry) ** 2
                           for entry, result in zip(entries, applied, strict=True))
    radius = ceil_sqrt(ceil_div(residual_squared, norm_squared))
    return Interval(center - radius, center + radius)


def model_spectral_error(matrix: list[list[int]], exponents: list[int]) -> tuple[dict, np.ndarray]:
    floating_matrix = np.array(matrix, dtype=float) / float(SCALE)
    values, vectors = eigh(floating_matrix)
    indices = np.argsort(np.abs(values))[::-1][:len(exponents)]
    assignments = [checked_eigenvalue(matrix, vectors[:, i]) for i in indices]
    ordered = sorted(assignments, key=lambda interval: interval.low)
    require(len(ordered) == len(exponents) and all(left.high < right.low
            for left, right in zip(ordered, ordered[1:])), "selected eigenvalue intervals overlap")

    omitted_squared = sum(entry ** 2 for row in matrix for entry in row)
    assigned_squared = maximum = 0
    for eigenvalue, exponent in zip(assignments, exponents, strict=True):
        singular = absolute(eigenvalue)
        ratio = amplitude(exponent)
        model = ratio * (ONE + ratio * ratio).reciprocal()
        omitted_squared -= singular.low ** 2
        assigned_squared += max(abs(singular.low - model.high), abs(singular.high - model.low)) ** 2
        maximum = max(maximum, model.high)
    require(omitted_squared >= 0, "assigned spectral mass exceeds Frobenius mass")
    return {
        "omitted_squared": omitted_squared,
        "assigned_squared": assigned_squared,
        "sorted_singular_error": ceil_sqrt(omitted_squared + assigned_squared),
        "maximum_model": maximum,
        "largest_eigenvalue_interval_width": max(x.high - x.low for x in assignments),
    }, values


def vector_distance(spectral_error: int, model_maximum: int) -> tuple[int, int, int]:
    maximum = spectral_error + model_maximum
    require(2 * maximum < SCALE, "the Schmidt angle derivative is unbounded on this enclosure")
    bound = Interval(maximum, maximum)
    derivative = (ONE - (bound * bound).times(4)).square_root().reciprocal()
    distance = SQRT_TWO * derivative * Interval(spectral_error, spectral_error)
    return distance.high, maximum, derivative.high


def certified_score(exponents: list[int], distance_upper: int) -> tuple[dict, Interval]:
    inverse_peak = ONE
    for exponent in exponents:
        term = ONE + amplitude(2 * exponent)
        inverse_peak = inverse_peak * term * term
    peak = inverse_peak.reciprocal()

    counts = rank_counts(exponents, ENERGY_CUTOFF)
    rank, mass, overlap, _ = run_statistics(counts)
    prefix_mass = peak * mass
    prefix_error = (prefix_mass.times(2) - SQRT_TWO * peak * overlap).nonnegative()
    tail = (ONE - prefix_mass).nonnegative()
    require(tail.high < ERROR_CAPS["omitted_geometric_probability_upper"],
            "omitted geometric probability exceeds the stated cap")
    distance = Interval(0, distance_upper)

    # The omitted overlap is nonnegative. The upper bound pays a linear
    # cost 2*tau in E for the unnormalized prefix, then transfers sqrt(E).
    model_upper = prefix_error + tail.times(2)
    physical_upper = model_upper.square_root() + distance.times(2)
    physical_lower = prefix_error.square_root() - (tail.square_root() + distance).times(2)
    physical_root = Interval(max(0, physical_lower.low), physical_upper.high)
    loss = (physical_root * physical_root).times(2).divide(243)
    return {
        "geometric_prefix_rank": rank,
        "geometric_prefix_mass": prefix_mass.report(),
        "geometric_prefix_duplication_defect": prefix_error.report(),
        "omitted_geometric_probability_upper": upper_text(tail.high),
        "geometric_duplication_defect_upper": upper_text(model_upper.high),
        "physical_duplication_defect": (physical_root * physical_root).report(),
    }, loss


def nominal_score(values: np.ndarray) -> dict:
    # Reuse the archive's unrounded mode-product enumeration for the sweep.
    squares = np.minimum(values ** 2, 0.25)
    occupations = 2 * squares / (1 + np.sqrt(1 - 4 * squares))
    probabilities = mode_probabilities(np.repeat(occupations, 2), NOMINAL_RELATIVE_THRESHOLD)
    defect = duplication_error(probabilities)
    return {
        "relative_probability_cutoff": NOMINAL_RELATIVE_THRESHOLD,
        "rank": len(probabilities),
        "retained_probability": float(probabilities.sum()),
        "duplication_defect": defect,
        "assisted_loss": 2 * defect / 243,
    }


def certify(power: int) -> dict:
    exponents = EXPONENTS[power]
    require(BASE == Interval.ratio(32767, 32768), "unexpected geometric amplitude ratio")
    require(len(exponents) == 24 and all(type(k) is int and k > 0 for k in exponents)
            and sorted(exponents) == exponents, "invalid geometric model exponents")
    size = 2 ** power
    compression = compression_bound(size).high
    matrix, rounding = enclosed_kernel(size)
    modes, values = model_spectral_error(matrix, exponents)
    singular_error = compression + rounding + modes["sorted_singular_error"]
    capped_errors = {
        "covariance_compression_singular_l2_upper": compression,
        "kernel_enclosure_frobenius_upper": rounding,
        "omitted_squared_singular_mass_upper": ceil_div(modes["omitted_squared"], SCALE),
        "retained_squared_singular_error_upper": ceil_div(modes["assigned_squared"], SCALE),
        "total_one_parity_singular_l2_upper": singular_error,
    }
    for name, upper in capped_errors.items():
        require(upper < ERROR_CAPS[name], f"{name} exceeds its stated cap at N = 2**{power}")
    distance, maximum, derivative = vector_distance(singular_error, modes["maximum_model"])
    distance_budget = SCORE_DISTANCE_BUDGETS[power]
    require(distance <= distance_budget, f"computed Schmidt distance exceeds the score budget at N = 2**{power}")
    scoring, loss = certified_score(exponents, distance_budget)
    below = loss.high < EXACT_FLOOR.low
    above = loss.low > EXACT_FLOOR.high
    require(below if power >= 19 else above, f"assisted loss at N = 2**{power} failed its strict comparison")
    reported_bound = Interval.ratio(REPORTED_LOSS_BOUNDS[power], 10 ** 6)
    require(loss.high < reported_bound.low if below else loss.low > reported_bound.high,
            f"assisted loss at N = 2**{power} failed its reported bound")
    gap = EXACT_FLOOR - loss if below else loss - EXACT_FLOOR
    return {
        "power": power,
        "chain_length": size,
        "compressed_covariance_dimension": len(matrix),
        "exponents": exponents,
        "nominal": nominal_score(values),
        "certified_error_bounds": {
            **{name: upper_text(value) for name, value in capped_errors.items()},
            "largest_selected_eigenvalue_interval_width_upper": upper_text(
                modes["largest_eigenvalue_interval_width"]),
            "rounded_matrix_singular_l2_upper": upper_text(modes["sorted_singular_error"]),
            "largest_physical_or_model_singular_value_upper": upper_text(maximum),
            "schmidt_angle_derivative_upper": upper_text(derivative),
            "sorted_schmidt_amplitude_l2_upper": upper_text(distance),
            "sorted_schmidt_amplitude_l2_budget_used": upper_text(distance_budget),
        },
        "certified_score": scoring,
        "certified_physical_assisted_loss": loss.report(),
        "comparison_to_flat_loss_floor": "below" if below else "above",
        "certified_absolute_gap": gap.report(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--powers", type=int, nargs="+", choices=tuple(EXPONENTS), default=list(EXPONENTS))
    args = parser.parse_args()
    if not __debug__:
        parser.error("The interval arithmetic checks require Python without -O")
    if len(set(args.powers)) != len(args.powers):
        parser.error("Each requested power must be unique")
    report = {
        "game": "original uniform 243-pair game",
        "shared_resource": "central-cut open XX ground state tensored with a shared Phi_9",
        "exact_flat_loss_floor": "(3-2*sqrt(2))/243",
        "certified_flat_loss_floor": EXACT_FLOOR.report(),
        "quadrature_order": ORDER,
        "geometric_amplitude_ratio": "32767/32768",
        "geometric_prefix_energy_cutoff": ENERGY_CUTOFF,
        "required_strict_error_caps": {name: upper_text(value) for name, value in ERROR_CAPS.items()},
        "rows": [certify(power) for power in args.powers],
    }
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
