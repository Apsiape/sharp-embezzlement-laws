"""Check the finite XX certificate's rounding bounds with exact rationals."""

from fractions import Fraction
import json
from math import factorial, isqrt
from pathlib import Path


SCALE = 10**45


def check_rounding_constants():
    unit = Fraction(1, 2**53)
    pi_upper = Fraction(22, 7)

    def gamma(steps):
        return steps * unit / (1 - steps * unit)

    jacobi_residual = Fraction("4e-14")
    jacobi_entry_perturbation = Fraction("4e-15")
    norm_slack = Fraction("2e-14")
    assert norm_slack**2 < (1 - norm_slack) * (jacobi_residual - jacobi_entry_perturbation)**2
    vector_error = Fraction("5.66e-11")
    gap = Fraction("1e-3") - jacobi_residual
    assert 2 * jacobi_residual**2 < ((vector_error - norm_slack) * gap)**2
    coordinate_error = 2 * vector_error + vector_error**2
    weight_error = Fraction("2e-7")
    assert ((coordinate_error + gamma(2) * (1 + norm_slack))
            / (Fraction("1e-3") - coordinate_error) < weight_error)
    node_error = Fraction("5e-14")
    assert (jacobi_residual + gamma(2) + gamma(3) / 2
            + (Fraction(1, 2) + jacobi_residual) * unit < node_error)

    pi_error = Fraction("1e-15") / 3
    argument_error = Fraction("3e-13")
    assert Fraction(1, 2**84 + 1) < unit
    assert Fraction(4, 3) * unit + unit * (1 + Fraction(4, 3) * unit) < 3 * unit
    assert ((1 + 4 * node_error / 3) * (1 + 3 * unit) * (1 + pi_error)
            * (1 + gamma(3)) - 1 < argument_error)
    node_difference_error = pi_upper * node_error / 2
    assert (node_difference_error + (pi_upper / 4 + node_difference_error)
            * ((1 + pi_error) * (1 + gamma(4)) - 1) < argument_error)

    for factorial_offset, argument_cap, argument_bound in (
        (1, Fraction("1.58"), pi_upper / 2 * (1 + argument_error)),
        (0, Fraction("0.8"), pi_upper / 4 + argument_error),
    ):
        assert argument_bound < argument_cap
        assert argument_cap**22 / factorial(22 + factorial_offset) < Fraction("1e-17")
        assert sum(argument_cap**(2 * j) / factorial(2 * j + factorial_offset) for j in range(11)) < 3
        ratio = argument_cap**2 / ((2 + factorial_offset) * (1 + factorial_offset))
        assert ratio * (1 + gamma(6)) < Fraction("0.42")
    assert factorial(21) < 2**66
    polynomial_error = Fraction("1e-13")
    trig_error = Fraction("1e-12")
    assert 3 * gamma(31) + Fraction("1e-17") < polynomial_error
    sinc_error = polynomial_error + pi_upper / 4 * argument_error
    assert ((1 + argument_error) * (1 + pi_upper / 2 * sinc_error) * (1 + unit)
            - 1 < trig_error)
    assert 2 * (polynomial_error + argument_error)**2 < trig_error**2
    reciprocal_error = Fraction("2e-12")
    assert (trig_error + unit) / (1 - trig_error) < reciprocal_error

    root_test_error = Fraction("1.1e-14")
    assert (1 + Fraction("1e-14")) * (1 + gamma(2)) - 1 < root_test_error
    weight_root_error = Fraction("3e-7")
    assert (1 - weight_root_error)**2 < (1 - weight_error)**2 * (1 - unit) * (1 - root_test_error)
    assert (1 + weight_error)**2 * (1 + unit) * (1 + root_test_error) < (1 + weight_root_error)**2
    entry_error = Fraction("1e-6")
    assert ((1 + reciprocal_error) * (1 + weight_root_error) * (1 + gamma(4))
            - 1 < entry_error)

    dimension = 1600
    trace_rounding = gamma(2 * dimension + 1)
    trace_error = Fraction("3e-12")
    assert 4 * trace_rounding / (1 - trace_rounding) < trace_error
    kernel_error = Fraction("2.1e-6")
    assert entry_error**2 * (4 + trace_error) < (kernel_error * (1 - entry_error))**2
    norm_lower = Fraction("0.99")
    underflow = Fraction(1, 2**1000)
    dot_rounding = gamma(dimension + 1)
    assert (40 * (2 * dimension * (1 + dot_rounding) + 1 + 1 / (1 - unit))
            / (norm_lower * 2**1075) < underflow)
    assert (1 - Fraction("1e-11") - underflow) / (1 + dot_rounding) > norm_lower**2
    assert 4 + trace_error < Fraction("2.01")**2
    assert (Fraction("4e-13") / (norm_lower * (1 - unit)) + Fraction("2.01") * dot_rounding
            + unit + underflow < Fraction("3e-12"))

    compression_error = Fraction("6e-6")
    assert Fraction(19, 4 * 3**18) + Fraction(4, 3 * 2**18) < compression_error
    assert Fraction(15, 4 * 3**24) + Fraction(4, 3 * 2**24) < Fraction("7.95e-8")
    singular_error = Fraction("4e-5") + compression_error + kernel_error
    assert singular_error == Fraction("4.81e-5")
    assert 16**2 * (1 - 4 * Fraction(".499")**2) > 1
    assert 2 * 16**2 < 23**2
    assert 23 * singular_error == Fraction("0.0011063") < Fraction(3, 2500)
    print("Binary64 rounding constants and Schmidt-vector transfer checked with exact rationals.")


def square_root_bounds(value):
    """Enclose a rational square root on the verifier's decimal grid."""
    scaled_square = value * SCALE**2
    lower = isqrt(scaled_square.numerator // scaled_square.denominator)
    upper = lower + (lower**2 < scaled_square)
    return Fraction(lower, SCALE), Fraction(upper, SCALE)


def check_score_transfers():
    sqrt_two_lower, sqrt_two_upper = square_root_bounds(Fraction(2))
    threshold_lower = Fraction(3, 2) - sqrt_two_upper
    threshold_upper = Fraction(3, 2) - sqrt_two_lower
    assert Fraction("0.0857864376") < threshold_lower <= threshold_upper < Fraction("0.0857864377")
    floor_lower = 2 * threshold_lower / 243
    floor_upper = 2 * threshold_upper / 243

    # The score verifier establishes this model bound. Use its rounded value
    # to check the physical bound without reading a prior score run.
    model_loss_upper = Fraction("0.000410720762")
    _, model_root_upper = square_root_bounds(model_loss_upper)
    physical_loss_upper = (model_root_upper + Fraction(3, 2500))**2
    assert physical_loss_upper < Fraction("0.000461") < floor_lower
    print("Chain-only transferred loss < 0.000461 < tracial floor.")

    report_path = Path(__file__).with_name("assisted_xx_certificate.json")
    rows = json.loads(report_path.read_text(encoding="utf-8"))["rows"]
    assert [row["power"] for row in rows] == [18, 19, 20]
    for row in rows:
        power = row["power"]
        assert row["chain_length"] == 2**power
        score = row["certified_score"]
        prefix = score["geometric_prefix_duplication_defect"]
        tail = Fraction(score["omitted_geometric_probability_upper"])
        distance = Fraction(row["certified_error_bounds"]["sorted_schmidt_amplitude_l2_budget_used"])
        prefix_root_lower, _ = square_root_bounds(Fraction(prefix["lower"]))
        _, tail_root_upper = square_root_bounds(tail)
        _, model_root_upper = square_root_bounds(Fraction(prefix["upper"]) + 2 * tail)
        defect_lower = max(0, prefix_root_lower - 2 * (distance + tail_root_upper))**2
        defect_upper = (model_root_upper + 2 * distance)**2

        defect = score["physical_duplication_defect"]
        recorded_lower, recorded_upper = Fraction(defect["lower"]), Fraction(defect["upper"])
        assert recorded_lower <= defect_lower <= defect_upper <= recorded_upper
        loss = row["certified_physical_assisted_loss"]
        loss_lower, loss_upper = Fraction(loss["lower"]), Fraction(loss["upper"])
        assert loss_lower <= 2 * recorded_lower / 243 <= 2 * recorded_upper / 243 <= loss_upper

        if power == 18:
            assert threshold_upper < recorded_lower and floor_upper < loss_lower
            print("Assisted 2^18 defect and loss exceed their thresholds.")
        else:
            stated_loss_upper = Fraction("0.000703" if power == 19 else "0.000663")
            assert recorded_upper < threshold_lower and loss_upper < stated_loss_upper < floor_lower
            print(f"Assisted 2^{power} defect and loss lie below their thresholds.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run the certificate without Python optimization (-O or PYTHONOPTIMIZE).")
    check_rounding_constants()
    check_score_transfers()

