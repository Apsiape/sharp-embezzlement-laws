"""Cross-check the explicit XX strategy with a shared Phi_9 at powers 19 and 20.

The primary verifier encloses the Jacobi roots and kernel entries directly.
This check instead bounds the fixed-point kernel perturbation from checked
Jacobi eigenpairs and uses a separate exact-rational mode-product model.
The pure-state loss bound applies to the model POVMs on the physical state.
The defect-based bound applies to the canonical POVMs in the physical
Schmidt basis. Both strategies use the same assisted shared resource.
"""

from fractions import Fraction
from math import factorial, isqrt, sqrt
import json

import numpy as np
from scipy.linalg import eigh, eigh_tridiagonal


SCALE = 10**28
RATIO_DENOMINATOR = 10**12
ORDER = 24
KEEP = 20
KERNEL_ERROR = Fraction(1, 10**8)
EIGEN_RESIDUAL = Fraction(1, 10**13)


def fixed(value):
    value = Fraction(value)
    return (value.numerator * SCALE) // value.denominator


def sqrt_enclosure(value):
    numerator = value.numerator * SCALE**2
    lower = isqrt(numerator // value.denominator)
    return Fraction(lower, SCALE), Fraction(lower + 1, SCALE)


def decimal_bound(value, *, upper):
    """Report a rational bound as a directed fixed-point decimal."""
    endpoint = fixed(value)
    if upper and Fraction(endpoint, SCALE) < value:
        endpoint += 1
    integer, fractional = divmod(endpoint, SCALE)
    return f"{integer}.{fractional:028d}"


def gauss_grid(m):
    """Propose and check Jacobi eigenpairs for the weighted quadrature.

    Each true Jacobi off-diagonal differs from its proposal by <1e-15.
    This costs <2e-15 in operator norm. The normalized spectral residual
    against the true Jacobi matrix is <3e-14, so the checked disjoint
    gaps >1e-3 identify all roots. Resolving the proposed vectors in the
    exact eigenbasis gives vector error <5e-11 after normalization and
    sign alignment. The checked first component then exceeds 1/32 in
    magnitude for the exact vector. The node relative error is <3.1e-14
    since L<=a, and the square-root weight relative error is <1.61e-9.
    Both bounds include fixed-point rounding.
    """
    nodes, weight_roots = [], []
    a = 1
    while a <= m:
        length = min(a, m - a + 1)
        if length <= ORDER:
            nodes.extend(i * SCALE for i in range(a, a + length))
            weight_roots.extend([SCALE] * length)
            a += length
            continue

        beta = []
        for j in range(1, ORDER):
            square = Fraction(j**2 * (length**2 - j**2), 4 * length**2 * (4*j**2 - 1))
            estimate = Fraction(sqrt(float(square)))
            error = Fraction(1, 10**15)
            assert 0 < estimate - error
            assert (estimate - error)**2 < square < (estimate + error)**2
            beta.append(estimate)

        values, vectors = eigh_tridiagonal(np.zeros(ORDER), np.array(beta, dtype=float))
        values = [Fraction(float(value)) for value in values]
        assert all(right - left > Fraction(1, 1000) for left, right in zip(values, values[1:]))
        length_root = isqrt(length * SCALE**2)
        for mu, vector in zip(values, vectors.T, strict=True):
            vector = [Fraction(float(value)) for value in vector]
            residual = [-mu * value for value in vector]
            for j, b in enumerate(beta):
                residual[j] += b * vector[j+1]
                residual[j+1] += b * vector[j]
            assert abs(sum(x*x for x in vector) - 1) < Fraction(2, 10**14)
            assert sum(x*x for x in residual) < Fraction(2, 10**14)**2
            assert vector[0]**2 > Fraction(1, 1000)
            nodes.append(fixed(Fraction(2*a + length - 1, 2) + length * mu))
            weight_roots.append(fixed(Fraction(length_root, SCALE) * abs(vector[0])))
        a += length
    return nodes, weight_roots


def atan_bounds(k):
    lower = sum(Fraction((-1)**j, (2*j+1)*k**(2*j+1)) for j in range(24))
    return lower, lower + Fraction(1, 49 * k**49)


def fixed_pi():
    low5, high5 = atan_bounds(5)
    low239, high239 = atan_bounds(239)
    low, high = 16*low5 - 4*high239, 16*high5 - 4*low239
    candidate = fixed(low)
    assert high < Fraction(candidate + 2, SCALE)
    return candidate


def polynomial(argument, coefficients):
    """Cosine or sinc of the fixed-point argument, rounded down at each step."""
    argument_square = argument * argument // SCALE
    value = coefficients[-1]
    for coefficient in reversed(coefficients[:-1]):
        value = coefficient + value * argument_square // SCALE
    return value


def rounded_kernel(size):
    """Evaluate the kernel at the fixed-point nodes and weight roots.

    The node bounds give a sine-argument relative error <5e-14 and a
    cosine-argument absolute error <5e-14, including the pi enclosure.
    The Horner rounding error is <2e-23 since the squared argument is
    below 2.5, with all fixed-point floors costing less than 1/SCALE.
    The verified Taylor remainders are <1e-17. The exact sinc and cosine
    exceed 2/pi and 1/sqrt(2), respectively. The resulting unweighted
    kernel relative error is <1e-12. The two square-root weight errors
    then give relative error <4e-9 in the matrix entries. The final
    fixed-point floors fit this allowance since the exact weights are
    >1/32 and the exact entries are >1/(32*(size+1)).
    The exact check ||rounded||_F<2 below therefore gives Frobenius
    perturbation <1e-8 from the exact quadrature kernel.
    """
    nodes, weight_roots = gauss_grid(size // 4)
    dim = len(nodes)
    matrix = [[0] * dim for _ in range(dim)]
    pi = fixed_pi()
    denominator = 2 * (size + 1)
    cosine = [fixed(Fraction((-1)**j, factorial(2*j))) for j in range(11)]
    sinc = [fixed(Fraction((-1)**j, factorial(2*j+1))) for j in range(11)]
    assert Fraction(158, 100)**22 / factorial(23) < Fraction(1, 10**17)
    assert Fraction(4, 5)**22 / factorial(22) < Fraction(1, 10**17)

    for i, x in enumerate(nodes):
        for j in range(i, dim):
            y = nodes[j]
            sum_angle = pi * (2*(x + y) - SCALE) // (SCALE * denominator)
            diff_angle = pi * (2*abs(x - y)) // (SCALE * denominator)
            assert 4 * SCALE // 10**6 < sum_angle < 158 * SCALE // 100
            assert 0 <= diff_angle < 4 * SCALE // 5
            sine_value = sum_angle * polynomial(sum_angle, sinc) // SCALE
            cosine_value = polynomial(diff_angle, cosine)
            cross = SCALE**2 // sine_value + SCALE**2 // cosine_value
            weight_product = weight_roots[i] * weight_roots[j] // SCALE
            matrix[i][j] = matrix[j][i] = cross * weight_product // (SCALE * denominator)
    return matrix


def mode_error(matrix):
    """Check exact squared residuals and bound matched and unmatched singular mass."""
    rounded = np.array(matrix, dtype=float) / float(SCALE)
    values, vectors = eigh(rounded)
    selection = np.argsort(np.abs(values))[::-1][:KEEP]
    norm_squared = Fraction(sum(entry**2 for row in matrix for entry in row), SCALE**2)
    assert norm_squared < 4
    assigned_lower = matched_upper = Fraction(0)
    ratios, proposals = [], []
    for i in selection:
        estimate = fixed(float(values[i]))
        vector = [fixed(float(x)) for x in vectors[:, i]]
        residual_squared = sum(
            (sum(entry * x for entry, x in zip(row, vector, strict=True)) - estimate * v)**2
            for row, v in zip(matrix, vector, strict=True)
        )
        vector_squared = sum(v*v for v in vector)
        assert vector_squared > 0
        assert residual_squared < (EIGEN_RESIDUAL * SCALE)**2 * vector_squared
        mu = Fraction(estimate, SCALE)
        proposals.append(mu)
        assert EIGEN_RESIDUAL < abs(mu) < Fraction(1, 2)
        s = float(abs(mu))
        r = 2*s / (1 + sqrt(1 - 4*s*s))
        a = round(RATIO_DENOMINATOR * r)
        assert 0 < a < RATIO_DENOMINATOR
        ratios.append(a)
        model_singular = Fraction(a * RATIO_DENOMINATOR, RATIO_DENOMINATOR**2 + a*a)
        matched_upper += (abs(model_singular - abs(mu)) + EIGEN_RESIDUAL)**2
        assigned_lower += (abs(mu) - EIGEN_RESIDUAL)**2

    proposals.sort()
    assert all(right - left > 2*EIGEN_RESIDUAL for left, right in zip(proposals, proposals[1:]))
    omitted = norm_squared - assigned_lower
    assert omitted >= 0
    return ratios, matched_upper, omitted


def prefix_bounds(ratios):
    """Sort exact coefficient numerators for the leading finite-model prefix."""
    Q = RATIO_DENOMINATOR
    peak_numerator = denominator = 1
    runs = [(1, 1)]
    for a in ratios:
        peak_numerator *= Q**2
        denominator *= Q**2 + a*a
        runs = [
            (coefficient * multiplier, count * multiplicity)
            for coefficient, count in runs
            for multiplier, multiplicity in ((Q**2, 1), (Q*a, 2), (a*a, 1))
            if coefficient * multiplier * 10**6 >= peak_numerator
        ]
    runs.sort(reverse=True)

    prefix_mass_numerator = rank = 0
    endpoints = []
    for a, count in runs:
        rank += count
        endpoints.append((rank, a))
        prefix_mass_numerator += count * a*a

    index = direct = dilated = overlap_numerator = 0
    while direct < len(endpoints):
        end, a = endpoints[direct]
        parent_end, b = endpoints[dilated]
        parent_end *= 2
        stop = min(end, parent_end)
        overlap_numerator += (stop - index) * a*b
        index = stop
        direct += stop == end
        dilated += stop == parent_end

    prefix_mass = Fraction(prefix_mass_numerator, denominator**2)
    bare_overlap = Fraction(overlap_numerator, denominator**2)
    sqrt2_low, _ = sqrt_enclosure(Fraction(2))
    assert prefix_mass < 1
    # The unit full model has E = 2 - sqrt(2) times its bare overlap.
    # The bare prefix overlap is a lower bound, since all terms are positive.
    return rank, len(runs), 1 - prefix_mass, 2 - sqrt2_low * bare_overlap


def check(size):
    if not __debug__:
        raise RuntimeError("Run the certificate with Python assertions enabled.")
    assert size in (2**19, 2**20)
    matrix = rounded_kernel(size)
    ratios, matched, omitted = mode_error(matrix)
    rank, runs, model_tail, model_error_upper = prefix_bounds(ratios)
    # Degree-23 dyadic polynomial compression bound from the main proof.
    compression = (
        Fraction(15, 4 * 3**ORDER)
        + Fraction(4, 3 * 2**ORDER)
    )
    assert compression < Fraction(8, 10**8)
    singular_error = sqrt_enclosure(matched + omitted)[1] + compression + KERNEL_ERROR
    peak_singular = max(Fraction(a * RATIO_DENOMINATOR, RATIO_DENOMINATOR**2 + a*a) for a in ratios)
    assert peak_singular + singular_error < Fraction(12, 25)
    _, sqrt2_high = sqrt_enclosure(Fraction(2))
    vector_error = sqrt2_high * Fraction(25, 7) * singular_error

    model_loss_upper = Fraction(2, 243) * model_error_upper
    generic_loss = (sqrt_enclosure(model_loss_upper)[1] + vector_error)**2
    specific_loss = Fraction(2, 243) * (sqrt_enclosure(model_error_upper)[1] + 2*vector_error)**2
    gamma_lower = (3 - 2*sqrt2_high) / 243
    stated_bound = Fraction(703 if size == 2**19 else 663, 10**6)
    assert specific_loss <= generic_loss < stated_bound < gamma_lower
    return {
        "size": size,
        "matrix_dimension": len(matrix),
        "mode_ratios": ratios,
        "ratio_denominator": RATIO_DENOMINATOR,
        "matched_error_squared_upper": decimal_bound(matched, upper=True),
        "omitted_singular_mass_upper": decimal_bound(omitted, upper=True),
        "one_parity_singular_error_upper": decimal_bound(singular_error, upper=True),
        "aligned_vector_error_upper": decimal_bound(vector_error, upper=True),
        "prefix_rank": rank,
        "prefix_runs": runs,
        "model_prefix_tail_upper": decimal_bound(model_tail, upper=True),
        "model_error_upper": decimal_bound(model_error_upper, upper=True),
        "model_loss_upper": decimal_bound(model_loss_upper, upper=True),
        "physical_loss_upper_generic_transfer": decimal_bound(generic_loss, upper=True),
        "physical_loss_upper_specific_transfer": decimal_bound(specific_loss, upper=True),
        "flat_loss_floor_lower": decimal_bound(gamma_lower, upper=False),
    }


if __name__ == "__main__":
    print(json.dumps({
        "game": "original uniform 243-pair game",
        "shared_resource": "central-cut open XX ground state tensored with Phi_9",
        "rows": [check(2**19), check(2**20)],
    }, indent=2))
