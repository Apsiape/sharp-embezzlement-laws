"""Certify the geometric mode model for the open XX chain at N = 2**84.

The NumPy/SciPy eigensolvers propose eigenpairs. The verifier uses residuals,
binary64 arithmetic bounds, and the dyadic quadrature estimate in the proof.
"""

from __future__ import annotations

from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, eigh_tridiagonal


SIZE = 2 ** 84
ORDER = 18
RATIO = Fraction(32767, 32768)
JACOBI_RESIDUAL = Fraction(4, 10 ** 14)
JACOBI_ENTRY_ERROR = Fraction(2, 10 ** 15)
CROSS_RESIDUAL = Fraction(3, 10 ** 12)
CROSS_TRACE_ERROR = Fraction(3, 10 ** 12)
COMPRESSION_ERROR = Fraction(6, 10 ** 6)
KERNEL_ERROR = Fraction(21, 10 ** 7)


def jacobi_pair_check(off_diagonal, values, orbitals):
    """Check each residual of the rounded Jacobi matrix in exact arithmetic."""
    off = [Fraction(float(entry)) for entry in off_diagonal]
    lambdas = [Fraction(float(value)) for value in values]
    assert min(b - a for a, b in zip(lambdas, lambdas[1:])) > Fraction(1, 1000)
    for value, vector in zip(lambdas, orbitals.T):
        v = [Fraction(float(entry)) for entry in vector]
        norm_squared = sum(entry ** 2 for entry in v)
        assert abs(norm_squared - 1) < Fraction(2, 10 ** 14)
        errors = [-value * entry for entry in v]
        for i, weight in enumerate(off):
            errors[i] += weight * v[i + 1]
            errors[i + 1] += weight * v[i]
        assert sum(error ** 2 for error in errors) < Fraction(2, 10 ** 14) ** 2
        assert v[0] ** 2 > Fraction(1, 1000)


def quadrature():
    nodes, weights = [], []
    first = 1
    count = SIZE // 4
    while first <= count:
        length = min(first, count - first + 1)
        if length <= ORDER:
            nodes.extend(float(i) for i in range(first, first + length))
            weights.extend([1.] * length)
        else:
            beta = []
            for k in range(1, ORDER):
                square = Fraction(k * k * (length ** 2 - k ** 2),
                                  4 * (4 * k ** 2 - 1) * length ** 2)
                candidate = math.sqrt(float(square))
                rounded = Fraction(candidate)
                assert rounded > JACOBI_ENTRY_ERROR
                assert (rounded - JACOBI_ENTRY_ERROR) ** 2 < square < (rounded + JACOBI_ENTRY_ERROR) ** 2
                beta.append(candidate)
            values, orbitals = eigh_tridiagonal(np.zeros(ORDER), np.array(beta))
            jacobi_pair_check(beta, values, orbitals)
            span = float(length)
            nodes.extend(float(first) + (span - 1) / 2 + span * values)
            weights.extend(span * orbitals[0, :] * orbitals[0, :])
        first += length
    return np.array(nodes), np.array(weights)


def arctan_enclosure(denominator):
    terms = [Fraction((-1) ** k, (2 * k + 1) * denominator ** (2 * k + 1)) for k in range(24)]
    lower = sum(terms)
    upper = lower + Fraction(1, 49 * denominator ** 49)
    return lower, upper


def verified_pi():
    lo5, hi5 = arctan_enclosure(5)
    lo239, hi239 = arctan_enclosure(239)
    lower, upper = 16 * lo5 - 4 * hi239, 16 * hi5 - 4 * lo239
    rounded = Fraction(math.pi)
    assert rounded - Fraction(1, 10 ** 15) < lower < upper < rounded + Fraction(1, 10 ** 15)
    return math.pi


def taylor_even(argument, first_factorial):
    """Degree 20 cosine or sinc Taylor polynomial, evaluated by Horner's rule."""
    square = argument * argument
    value = np.full_like(square, 1 / math.factorial(20 + first_factorial))
    for k in reversed(range(10)):
        value = value * square + (-1) ** k / math.factorial(2 * k + first_factorial)
    return value


def rounded_cross_matrix():
    nodes, weights = quadrature()
    denominator = float(2 * (SIZE + 1))
    pi = verified_pi()
    twice = 2 * nodes
    sum_argument = pi * (twice[:, None] + twice[None, :] - 1) / denominator
    difference_argument = pi * np.abs(twice[:, None] - twice[None, :]) / denominator
    cross = 1 / (sum_argument * taylor_even(sum_argument, 1))
    cross += 1 / taylor_even(difference_argument, 0)
    weight_products = weights[:, None] * weights[None, :]
    weight_roots = np.sqrt(weight_products)
    assert np.min(weight_roots) > 0
    assert np.max(np.abs(weight_roots * weight_roots / weight_products - 1)) < 1e-14
    cross *= weight_roots / denominator
    assert np.array_equal(cross, cross.T)
    assert len(cross) < 1600
    return cross


def decimal_power_interval(exponent):
    endpoints = []
    for rounding in (ROUND_FLOOR, ROUND_CEILING):
        with localcontext() as context:
            context.prec = 40
            context.rounding = rounding
            base = Decimal(RATIO.numerator) / Decimal(RATIO.denominator)
            result = Decimal(1)
            remaining = exponent
            while remaining:
                if remaining % 2:
                    result *= base
                base *= base
                remaining //= 2
            endpoints.append(Fraction(result))
    return endpoints


def model_singular_interval(exponent):
    lower, upper = decimal_power_interval(exponent)
    return lower / (1 + lower ** 2), upper / (1 + upper ** 2)


def certify_rounded_cross(matrix, exponents):
    """Bound the distance between rounded-kernel singular values and the model."""
    trace_estimate = float(np.sum(np.sum(matrix * matrix, axis=1)))
    assert trace_estimate < 4
    trace_upper = Fraction(trace_estimate) + CROSS_TRACE_ERROR
    eigenvalues, eigenvectors = eigh(matrix)
    selected = np.argsort(np.abs(eigenvalues))[::-1][:len(exponents)]
    signed_values = sorted(Fraction(float(eigenvalues[i])) for i in selected)
    assert min(b - a for a, b in zip(signed_values, signed_values[1:])) > Fraction(1, 10 ** 10)

    assigned_mass_lower = Fraction(0)
    selected_squared_error = Fraction(0)
    max_model_value = Fraction(0)
    for i, exponent in zip(selected, exponents):
        value, vector = eigenvalues[i], eigenvectors[:, i]
        assert abs(value) < 1
        assert abs(float(np.sum(vector * vector)) - 1) < 1e-11
        residual = np.sum(matrix * vector[None, :], axis=1) - value * vector
        assert np.max(np.abs(residual)) < 1e-14
        absolute = abs(Fraction(float(value)))
        assert absolute > CROSS_RESIDUAL
        assigned_mass_lower += (absolute - CROSS_RESIDUAL) ** 2
        lower, upper = model_singular_interval(exponent)
        max_model_value = max(max_model_value, upper)
        selected_squared_error += (max(abs(absolute - lower), abs(absolute - upper)) + CROSS_RESIDUAL) ** 2

    omitted_squared_error = trace_upper - assigned_mass_lower
    assert omitted_squared_error >= 0
    return omitted_squared_error, selected_squared_error, max_model_value


def main():
    exponents_path = Path(__file__).with_name("model84_exponents.json")
    exponents = json.loads(exponents_path.read_text())
    assert len(exponents) == 81 and all(isinstance(value, int) and value > 0 for value in exponents)
    matrix = rounded_cross_matrix()
    assert len(matrix) == 1418
    omitted, quantized, maximum = certify_rounded_cross(matrix, exponents)
    rounded_singular_error = Fraction(1, 25000)
    assert omitted + quantized < rounded_singular_error ** 2
    singular_error = rounded_singular_error + COMPRESSION_ERROR + KERNEL_ERROR
    assert maximum + singular_error < Fraction(499, 1000)
    # The doubled Gaussian state is a product of two copies of the mode list.
    # On [0, .499], the map from s to its Schmidt angle has derivative < 16.
    vector_error = Fraction(23) * singular_error
    assert vector_error < Fraction(12, 10000)
    print(json.dumps(dict(size=SIZE, nodes=len(matrix), modes=len(exponents),
                          rounded_singular_error=float(rounded_singular_error),
                          omitted_squared_error=float(omitted),
                          quantization_squared_error=float(quantized),
                          singular_error=float(singular_error),
                          vector_error=float(vector_error)), indent=2))


if __name__ == "__main__":
    main()
