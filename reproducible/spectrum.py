"""Spectral formulas, ordered rank runs, and truncation estimates."""

from itertools import accumulate
import math

import numpy as np

from game import CLASSICAL_LOSS, strategy_loss


LOG2 = math.log(2)
FLAT_ERROR = 2 - math.sqrt(2)
ASSISTED_COEFFICIENT = 2 / 243
GLOBAL_BLOCK_BOUND = (
    math.sqrt(ASSISTED_COEFFICIENT)
    + (1 + 2 * math.sqrt(ASSISTED_COEFFICIENT)) * math.sqrt((8 / 3) * (6 + 4 * math.sqrt(2)))
) ** 2


def error(probabilities, dilation=2):
    """Squared dilation defect of a decreasing list, homogeneous in its mass."""
    coefficients = np.sqrt(probabilities)
    expanded = np.repeat(coefficients, dilation) / math.sqrt(dilation)
    expanded[:len(coefficients)] -= coefficients
    return float(np.dot(expanded, expanded))


def compressed_error(log_weights, multiplicities, dilation=2):
    """The same defect for consecutive equal-weight runs, using exact rank counts."""
    runs = [(float(weight), int(count)) for weight, count in zip(log_weights, multiplicities, strict=True)
            if count]
    if not runs:
        return 0.0
    log_weights, counts = zip(*runs)
    ends = list(accumulate(counts))
    log_dilation = math.log(dilation)
    left = 0
    direct, dilated = 0, 0
    terms = []
    while dilated < len(ends):
        direct_end = ends[direct] if direct < len(ends) else dilation * ends[-1]
        dilated_end = dilation * ends[dilated]
        right = min(direct_end, dilated_end)
        log_direct = log_weights[direct] if direct < len(ends) else -math.inf
        log_dilated = log_weights[dilated] - log_dilation
        larger, smaller = max(log_direct, log_dilated), min(log_direct, log_dilated)
        terms.append(math.exp(math.log(right - left) + larger)
                     * math.expm1((smaller - larger) / 2) ** 2)
        left = right
        direct += right == direct_end
        dilated += right == dilated_end
    return math.fsum(terms)


def block_average(probabilities, size=9):
    padded = np.pad(probabilities, (0, (-len(probabilities)) % size))
    averaged = np.repeat(padded.reshape(-1, size).mean(axis=1), size)
    difference = np.sqrt(padded) - np.sqrt(averaged)
    return averaged, float(np.dot(difference, difference))


def sqrt_tail_bound(tau):
    """Bound on |sqrt(E) - sqrt(E of normalized retained probabilities)|."""
    # Use 1 - sqrt(1 - tau) = tau / (1 + sqrt(1 - tau)) to avoid cancellation.
    return 2 * math.sqrt(2 * tau / (1 + math.sqrt(1 - tau)))


def tail_interval(normalized_error, tau, upper_limit=FLAT_ERROR):
    center = math.sqrt(normalized_error)
    radius = sqrt_tail_bound(tau)
    return max(0., center - radius) ** 2, min(upper_limit, (center + radius) ** 2)


def mode_probabilities(modes, relative_threshold=1e-11, cap=12_000_000):
    """Enumerate Bernoulli products whose probability exceeds a relative threshold."""
    modes = np.clip(modes, 0., 0.5)
    peak = math.exp(float(np.log1p(-modes).sum()))
    ratios = np.sort(modes / (1 - modes))[::-1]
    products = np.ones(1)
    for ratio in ratios[ratios > relative_threshold]:
        occupied = products * ratio
        products = np.concatenate([products, occupied[occupied > relative_threshold]])
        if len(products) > cap:
            raise ValueError("The product enumeration exceeded the cap. Increase the relative threshold.")
    return np.sort(products)[::-1] * peak


def rank_moments(probabilities):
    ranks = np.arange(len(probabilities) + 1, dtype=float)
    logs = np.zeros_like(ranks)
    np.log(ranks, out=logs, where=ranks > 0)
    first = np.dot(probabilities, np.diff(ranks * (logs - 1)))
    second = np.dot(probabilities, np.diff(ranks * (logs ** 2 - 2 * logs + 2)))
    return float((second - first ** 2) / LOG2 ** 2)


def spectral_row(probabilities, missing_mode_mass_bound=None):
    """Diagnostics for retained probabilities, including the two explicit strategies."""
    probabilities = np.asarray(probabilities, dtype=float)
    probabilities = probabilities[probabilities > 0]
    logs = np.log(probabilities)
    mass = float(probabilities.sum())
    normalized = probabilities / mass
    lost_mass = max(0., 1 - mass)
    e = error(probabilities)
    e_normalized = e / mass
    e_low, e_high = tail_interval(e_normalized, lost_mass)
    averaged, distance_squared = block_average(normalized)
    e_averaged = error(averaged)
    mean_log = np.dot(probabilities, logs)
    row = dict(
        S=float(-mean_log), b=-float(logs[0]), E=e, varX=rank_moments(probabilities),
        varK=float(np.dot(probabilities, (logs - mean_log) ** 2)), mass=mass, n=len(probabilities),
        E_normalized=e_normalized, E_lower=e_low, E_upper=e_high,
        tail_mass=lost_mass, sqrt_E_tail_radius=sqrt_tail_bound(lost_mass),
        d9_squared=distance_squared, E_B9=e_averaged,
        loss_B9=(math.sqrt(distance_squared) + math.sqrt(ASSISTED_COEFFICIENT * e_averaged)) ** 2,
        loss_explicit=strategy_loss(normalized),
    )
    if missing_mode_mass_bound is not None:
        total_tau = 1 - mass * (1 - missing_mode_mass_bound)
        total_tau = min(1., max(0., total_tau))
        e_mode_low, e_mode_high = tail_interval(e_normalized, total_tau)
        row.update(
            missing_mode_mass_bound=missing_mode_mass_bound,
            tail_mass_with_modes_bound=total_tau,
            E_lower_with_modes=e_mode_low, E_upper_with_modes=e_mode_high,
            loss_upper_prefix=(1 - total_tau) * min(CLASSICAL_LOSS, row["loss_explicit"])
            + total_tau * CLASSICAL_LOSS,
        )
    return row
