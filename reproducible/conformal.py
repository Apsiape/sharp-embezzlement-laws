"""Conformal tower spectra, the smooth CL density, and repeated pair spectra."""

from functools import lru_cache
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import i0e

from spectrum import compressed_error, tail_interval


def tower_counts(copies, max_energy):
    """Coefficients of the odd-level product, with exact integer degeneracies."""
    counts = [1] + [0] * max_energy
    for level in range(1, max_energy + 1, 2):
        for _ in range(copies):
            for energy in range(max_energy, level - 1, -1):
                counts[energy] += counts[energy - level]
    return counts


def tower_partition(epsilon, copies):
    """Finite sums for log Z and entropy with analytic positive tail bounds."""
    levels = range(1, math.ceil(45 / epsilon) + 1, 2)
    energies = epsilon * np.array(levels)
    weights = np.exp(-energies)
    log_z = float(copies * np.log1p(weights).sum())
    entropy = log_z + float(copies * np.dot(energies, weights / (1 + weights)))
    next_energy = (levels[-1] + 2) * epsilon
    next_weight = math.exp(-next_energy)
    ratio = math.exp(-2 * epsilon)
    gap = -math.expm1(-2 * epsilon)
    z_tail = copies * next_weight / gap
    entropy_tail = z_tail * (1 + next_energy + 2 * epsilon * ratio / gap)
    return log_z, entropy, z_tail, entropy_tail


def tower_parameters(epsilon, copies, tail_target):
    log_z, entropy, z_tail, entropy_tail = tower_partition(epsilon, copies)
    log_z_half, _, half_tail, _ = tower_partition(epsilon / 2, copies)
    # Markov bound with tilting parameter epsilon / 2.
    prefactor = log_z_half + half_tail - log_z
    cutoff = math.ceil(2 * (prefactor - math.log(tail_target)) / epsilon)
    tau_bound = math.exp(prefactor - epsilon * (cutoff + 1) / 2)
    return cutoff, log_z, entropy, tau_bound, z_tail, entropy_tail


def tower_rows(epsilons, copies, dilation=2, tail_target=1e-15):
    parameters = [tower_parameters(epsilon, copies, tail_target) for epsilon in epsilons]
    counts = tower_counts(copies, max(row[0] for row in parameters))
    flat_bound = 2 - 2 / math.sqrt(dilation)
    for epsilon, (cutoff, log_z, entropy, tau_bound, z_tail, entropy_tail) in zip(epsilons, parameters):
        log_weights = [-log_z - epsilon * n for n in range(cutoff + 1)]
        kept = counts[:cutoff + 1]
        mass = math.fsum(math.exp(math.log(count) + weight)
                         for count, weight in zip(kept, log_weights) if count)
        value = compressed_error(log_weights, kept, dilation) / mass
        low, high = tail_interval(value, tau_bound, flat_bound)
        yield dict(copies=copies, eps=epsilon, b=log_z, S=entropy, E=value, mass=mass,
                   dilation=dilation, E_lower=low, E_upper=high, energy_cutoff=cutoff,
                   probability_tail_bound=tau_bound, log_z_tail_bound=z_tail,
                   entropy_tail_bound=entropy_tail)


def log_i0(z):
    return z + math.log(i0e(z))


@lru_cache(maxsize=16384)
def inverse_log_i0(log_rank):
    if log_rank <= 0:
        return 0.
    upper = log_rank + 4 * math.sqrt(log_rank + 1)
    while log_i0(upper) < log_rank:
        upper *= 2
    return brentq(lambda z: log_i0(z) - log_rank, 0, upper, xtol=2e-13)


def smooth_row(b, dilation=2, tail_target=1e-15):
    """Integrate the squared log-rank translation defect of the smooth CL density."""
    h = math.log(dilation)
    if b <= 0 or h <= 0:
        raise ValueError("Use b > 0 and dilation > 1.")

    def log_amplitude(y):
        z = inverse_log_i0(y)
        return (y - b - z ** 2 / (4 * b)) / 2

    def difference_squared(y):
        a, c = log_amplitude(y), log_amplitude(y - h)
        larger, smaller = max(a, c), min(a, c)
        return math.exp(2 * larger) * math.expm1(smaller - larger) ** 2

    # The CL law for T=z²/(4b) is compound Poisson. Its Chernoff bound above
    # z_cut is exp(-(z_cut-2b)²/(4b)), and the second integrand is shifted by h.
    z_cut = 2 * b + 2 * math.sqrt(-b * math.log(tail_target))
    y_cut = log_i0(z_cut)
    upper = y_cut + h
    intervals = sorted({0., h, max(h, log_i0(2 * b)), upper})
    value = math.exp(-b) * math.expm1(-h / 2) ** 2
    estimated_error = 0.
    for left, right in zip(intervals, intervals[1:]):
        part, quad_error = quad(difference_squared, left, right, epsabs=1e-13, epsrel=1e-12)
        value += part
        estimated_error += quad_error

    density = lambda y: math.exp(2 * log_amplitude(y))
    mass, mass_error = quad(density, 0., y_cut, epsabs=1e-13, epsrel=1e-12)
    return dict(b=b, S=2 * b, E=value, dilation=dilation, mass=mass + math.exp(-b),
                E_right_tail_bound=2 * tail_target, quad_error_estimate=estimated_error,
                mass_quad_error_estimate=mass_error)


def binomial_counts(number):
    counts = [1]
    for index in range(number):
        counts.append(counts[-1] * (number - index) // (index + 1))
    return counts


def pair_row(larger_weight, number, dilation=2):
    p = larger_weight
    lp, lq = math.log(p), math.log1p(-p)
    weights = [(number - k) * lp + k * lq for k in range(number + 1)]
    value = compressed_error(weights, binomial_counts(number), dilation)
    return dict(lam=p, n=number, E=value, S=-number * (p * lp + (1 - p) * lq), dilation=dilation)
