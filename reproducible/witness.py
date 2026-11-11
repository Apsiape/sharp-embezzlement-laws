"""Verify the compact and dyadic sine witnesses for the 243-pair game."""

import argparse
from fractions import Fraction
import json
import math

from certificate import Quadratic
from game import QUESTIONS_PER_GAME, TRACIAL_LOSS_FLOOR, strategy_loss
from spectrum import ASSISTED_COEFFICIENT, compressed_error


def compact_witness():
    shells = (15, 11, 10, 8, 6, 3)
    norm_squared = sum(value**2 for value in shells)
    assert norm_squared == 555

    amplitudes = []
    for shell, value in enumerate(shells):
        exponent = max(0, shell - 1)
        rational_scale = Fraction(value, 2 ** (exponent // 2))
        coefficient = (
            Quadratic(b=rational_scale / 2) if exponent % 2 else Quadratic(rational_scale)
        )
        amplitudes.extend([coefficient] * (2 ** exponent))

    # Repeat and rescale the amplitudes for the shared Phi_9, keeping
    # the common normalization sqrt(555) outside the exact calculation.
    repeated = [value * Quadratic(Fraction(1, 3)) for value in amplitudes for _ in range(9)]
    weights = [value * value for value in repeated]
    assert all(weight.b == 0 for weight in weights)
    assert sum(weight.a for weight in weights) == norm_squared
    assert all(left.a >= right.a for left, right in zip(weights, weights[1:]))

    inverse_sqrt_two = Quadratic(b=Fraction(1, 2))
    padded = repeated + [Quadratic()] * len(repeated)
    differences = [
        value - inverse_sqrt_two * repeated[index // 2]
        for index, value in enumerate(padded)
    ]
    exact_error = sum((value * value for value in differences), Quadratic())
    exact_error = exact_error * Quadratic(Fraction(1, norm_squared))
    assert exact_error == Quadratic(Fraction(598, 555), Fraction(-390, 555))

    exact_loss = exact_error * Quadratic(Fraction(2, QUESTIONS_PER_GAME))
    exact_floor = Quadratic(Fraction(3, QUESTIONS_PER_GAME), Fraction(-2, QUESTIONS_PER_GAME))
    exact_margin = exact_floor - exact_loss
    denominator = norm_squared * QUESTIONS_PER_GAME
    assert exact_margin == Quadratic(Fraction(469, denominator), Fraction(-330, denominator))
    assert exact_margin.positive()

    error = float(exact_error.a) + float(exact_error.b) * math.sqrt(2)
    loss = float(exact_loss.a) + float(exact_loss.b) * math.sqrt(2)
    probabilities = [float(weight.a / norm_squared) for weight in weights]
    assert math.isclose(strategy_loss(probabilities), loss, rel_tol=0, abs_tol=1e-15)
    rank = len(repeated)
    return dict(
        shell_coefficients=list(shells), shell_norm_squared=norm_squared,
        base_schmidt_rank=len(amplitudes), schmidt_rank=rank, qubits_per_party=(rank - 1).bit_length(),
        E=error, E_exact={"rational": str(exact_error.a), "sqrt2_coefficient": str(exact_error.b)},
        explicit_loss=loss, tracial_loss_floor=TRACIAL_LOSS_FLOOR,
        certified_margin=float(exact_margin.a) + float(exact_margin.b) * math.sqrt(2),
    )


def sine_witness(shells=9):
    theta = math.pi / (shells + 3)
    x = [math.sin((k + 2) * theta) for k in range(1, shells + 1)]
    x.insert(0, math.sqrt(2) * x[0])
    counts = [1] + [2 ** (k - 1) for k in range(1, shells + 1)]
    norm_squared = math.fsum(value ** 2 for value in x)
    p = [value ** 2 / (norm_squared * count) for value, count in zip(x, counts, strict=True)]
    assert all(a >= b for a, b in zip(p, p[1:]))
    shell_formula = ((math.sqrt(2) - 1) ** 2 * x[1] ** 2
                     + math.fsum((right - left) ** 2 for left, right in zip(x[1:], x[2:] + [0.]))) / norm_squared
    ordinary = compressed_error([math.log(value) for value in p], counts)
    repeated = compressed_error([math.log(value / 9) for value in p], [9 * count for count in counts])
    assert math.isclose(ordinary, shell_formula, abs_tol=1e-14)
    assert math.isclose(repeated, ordinary, abs_tol=1e-14)
    upper = ASSISTED_COEFFICIENT * theta ** 2
    return dict(
        shells=shells, schmidt_rank=9 * 2 ** shells, qubits_per_party=(9 * 2 ** shells - 1).bit_length(),
        E=ordinary, E_shell_sum=shell_formula, E_repeated=repeated,
        explicit_loss=ASSISTED_COEFFICIENT * repeated,
        proved_loss_upper=upper,
        tracial_loss_floor=TRACIAL_LOSS_FLOOR,
        certified_margin=TRACIAL_LOSS_FLOOR - upper,
        loss_upper_after_vector_error=(math.sqrt(upper) + 0.002) ** 2,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compact", action="store_true", help="Verify the nine-qubit witness.")
    args = parser.parse_args()
    print(json.dumps(compact_witness() if args.compact else sine_witness(), indent=2))
