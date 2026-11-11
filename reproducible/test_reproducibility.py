"""Numerical regression and independent checks for the manuscript supplement."""

from itertools import product
import json
import math
from pathlib import Path
import unittest

import numpy as np

from certificate import verify_certificate
from chains import chain_row, dense_modes, ising_modes
from conformal import pair_row, smooth_row, tower_counts, tower_rows
from game import (CLASSICAL_LOSS, QUESTIONS_PER_GAME, classical_strategy_lost_pairs,
                  direct_strategy_loss, fms_loss_kernel, measurements, strategy_loss)
from spectrum import (GLOBAL_BLOCK_BOUND, LOG2, block_average, compressed_error, error,
                      mode_probabilities, rank_moments, tail_interval)
from witness import compact_witness, sine_witness


ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT / "original" / "supplement" / "data"
GENERATED = ROOT / "generated"


def spin_hamiltonian(model, size, delta=0.):
    one = np.eye(2)
    x = np.array([[0., 1.], [1., 0.]])
    y = np.array([[0., -1j], [1j, 0.]])
    z = np.diag([1., -1.])

    def place(observable, site):
        result = np.eye(1)
        for i in range(size):
            result = np.kron(result, observable if i == site else one)
        return result

    xs, ys, zs = [[place(operator, j) for j in range(size)] for operator in (x, y, z)]
    if model == "ising":
        return -sum(xs[j] @ xs[j + 1] for j in range(size - 1)) - sum(zs)
    xx = -sum((xs[j] @ xs[j + 1] + ys[j] @ ys[j + 1]).real / 2 for j in range(size - 1))
    return xx + delta * sum((-1) ** j * (np.eye(2 ** size) - zs[j]) / 2 for j in range(size))


def schmidt_from_spin_hamiltonian(hamiltonian):
    _, states = np.linalg.eigh(hamiltonian)
    dimension = math.isqrt(len(states))
    return np.linalg.svd(states[:, 0].reshape(dimension, dimension), compute_uv=False) ** 2


class Reproducibility(unittest.TestCase):
    def test_dilation_and_compressed_ranks(self):
        rng = np.random.default_rng(12)
        for dilation in (2, 3, 9):
            probabilities = np.repeat(np.sort(rng.random(7))[::-1], [2, 4, 1, 6, 8, 5, 2])
            probabilities /= probabilities.sum()
            coefficients = np.sqrt(probabilities)
            indices = np.arange(len(coefficients)) // dilation
            overlap = np.dot(coefficients, coefficients[indices]) / math.sqrt(dilation)
            self.assertAlmostEqual(error(probabilities, dilation), 2 - 2 * overlap, places=14)
            distinct, counts = np.unique(probabilities, return_counts=True)
            self.assertAlmostEqual(compressed_error(np.log(distinct[::-1]), counts[::-1], dilation),
                                   error(probabilities, dilation), places=14)
            replicated = np.repeat(probabilities, 7) / 7
            self.assertAlmostEqual(error(replicated, dilation), error(probabilities, dilation), places=14)

    def test_bounds_and_tail(self):
        rng = np.random.default_rng(2026)
        samples = [np.full(n, 1 / n) for n in (1, 4, 9, 19)]
        samples += [1 / np.arange(1, n + 1) for n in (20, 80, 250, 1000)]
        samples += [np.sort(rng.random(n) ** rng.uniform(0.1, 12))[::-1] for n in range(2, 137, 3)]
        for p in samples:
            p = p / p.sum()
            e = error(p)
            s = -np.dot(p, np.log(p))
            var = rank_moments(p)
            self.assertGreaterEqual(e, p[0] / (6 + 4 * math.sqrt(2)) - 1e-14)
            self.assertGreaterEqual(math.sqrt(e), LOG2 / (2 * (s + 1 + 3 * LOG2)))
            self.assertGreaterEqual(e * (4 * var + 1), 1)
            q, d_squared = block_average(p)
            self.assertLessEqual(d_squared, 8 * p[0] / 3 + 1e-14)
            adaptive = (math.sqrt(d_squared) + math.sqrt(2 * error(q) / QUESTIONS_PER_GAME)) ** 2
            self.assertLessEqual(adaptive, GLOBAL_BLOCK_BOUND * e + 1e-14)
            self.assertLessEqual(strategy_loss(p), (7 / 3) * e + 1e-14)
            kept = p[:max(1, 3 * len(p) // 4)]
            low, high = tail_interval(error(kept / kept.sum()), max(0., 1 - kept.sum()))
            self.assertLessEqual(low, e + 1e-14)
            self.assertGreaterEqual(high, e - 1e-14)
        _, distance = block_average(np.full(4, 0.25))
        self.assertAlmostEqual(distance, 2 / 3, places=14)

    def test_spin_exact_diagonalization(self):
        for size, (model, delta) in product((4, 8), (("ising", 0.), ("xx", 0.), ("xx", 0.17))):
            with self.subTest(size=size, model=model, delta=delta):
                expected = schmidt_from_spin_hamiltonian(spin_hamiltonian(model, size, delta))
                modes = ising_modes(size).occupation if model == "ising" else dense_modes(size, delta)
                spectrum = mode_probabilities(modes, relative_threshold=1e-18)
                spectrum = np.pad(spectrum, (0, len(expected) - len(spectrum)))
                np.testing.assert_allclose(expected, spectrum, atol=4e-12, rtol=0)

    def test_dense_fast_and_doubled_chains(self):
        size = 512
        dense = mode_probabilities(dense_modes(size))
        direct = chain_row("xx", size, count=32)
        doubled = chain_row("xx-doubled", size, count=32)
        self.assertAlmostEqual(error(dense), direct["E"], places=12)
        self.assertAlmostEqual(doubled["E"], direct["E"], places=12)
        self.assertAlmostEqual(doubled["S"], direct["S"], places=12)

    def test_fms_measurements_and_loss(self):
        for side in measurements():
            for effects in side.values():
                np.testing.assert_allclose(sum(effects.values()), np.eye(5), atol=2e-15)
                for effect in effects.values():
                    np.testing.assert_allclose(effect @ effect, effect, atol=2e-15)
        expected = np.full((5, 5), -1.5)
        expected[0, 1:] = expected[1:, 0] = -math.sqrt(2) / 4
        expected[1, 2] = expected[2, 1] = expected[3, 4] = expected[4, 3] = -2.375
        np.fill_diagonal(expected, [2, 5.625, 5.625, 5.625, 5.625])
        np.testing.assert_allclose(fms_loss_kernel(), expected, atol=1e-14)
        rng = np.random.default_rng(17)
        for rank in (1, 2, 4, 8, 9, 10, 17, 24, 47):
            p = np.sort(rng.random(rank))[::-1]
            p /= p.sum()
            self.assertAlmostEqual(direct_strategy_loss(p), strategy_loss(p), places=14)
        self.assertEqual(len(classical_strategy_lost_pairs()) / QUESTIONS_PER_GAME, CLASSICAL_LOSS)
        self.assertEqual(len(verify_certificate()), 19)

    def test_theory_and_pairs(self):
        counts = tower_counts(1, 6)
        self.assertEqual(counts, [1, 1, 0, 1, 1, 1, 1])
        for p in (2 / 3, 0.75, 0.8):
            n = 10
            expanded = np.repeat([p ** (n - k) * (1 - p) ** k for k in range(n + 1)],
                                 [math.comb(n, k) for k in range(n + 1)])
            self.assertAlmostEqual(pair_row(p, n)["E"], error(expanded), places=13)
        for row in tower_rows([0.5, 0.05], 2):
            self.assertAlmostEqual(row["mass"], 1, places=12)
            self.assertLess(row["probability_tail_bound"], 1e-15)
        cl = smooth_row(64)
        self.assertAlmostEqual(cl["mass"], 1, places=12)
        self.assertAlmostEqual(cl["E"], 0.0009418725133298889, places=10)

    def test_dyadic_sine_witness(self):
        witness = sine_witness()
        self.assertEqual(witness["shells"], 9)
        self.assertEqual(witness["schmidt_rank"], 4608)
        self.assertEqual(witness["qubits_per_party"], 13)
        self.assertAlmostEqual(witness["explicit_loss"], 0.000402257840, places=12)
        self.assertAlmostEqual(witness["proved_loss_upper"], 0.000564106333, places=12)
        self.assertAlmostEqual(witness["tracial_loss_floor"], 0.000706061215, places=12)
        self.assertLess(witness["explicit_loss"], witness["proved_loss_upper"])
        self.assertLess(witness["loss_upper_after_vector_error"], witness["tracial_loss_floor"])

    def test_compact_witness(self):
        witness = compact_witness()
        self.assertEqual(witness["shell_coefficients"], [15, 11, 10, 8, 6, 3])
        self.assertEqual(witness["shell_norm_squared"], 555)
        self.assertEqual(witness["base_schmidt_rank"], 32)
        self.assertEqual(witness["schmidt_rank"], 288)
        self.assertEqual(witness["qubits_per_party"], 9)
        self.assertEqual(witness["E_exact"], {"rational": "598/555", "sqrt2_coefficient": "-26/37"})
        self.assertAlmostEqual(witness["explicit_loss"], 0.000688936502, places=12)
        self.assertLess(witness["explicit_loss"], witness["tracial_loss_floor"])

    def test_supplied_critical_points(self):
        critical = json.loads((GENERATED / "critical.json").read_text())
        for model, filename in (("ising", "data_ising_3_17.json"), ("xx-doubled", "data_xx_sq.json")):
            archive = {row["size"]: row for row in json.loads((ARCHIVE / filename).read_text())}
            sample = [row for row in critical if row["model"] == model]
            self.assertEqual({row["size"] for row in sample}, set(archive))
            for row in sample:
                original = archive[row["size"]]
                for field in ("S", "b", "E", "varX", "varK", "mass"):
                    self.assertAlmostEqual(row[field], original[field], places=9)


if __name__ == "__main__":
    unittest.main(verbosity=2)
