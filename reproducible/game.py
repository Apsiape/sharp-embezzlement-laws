"""Loss of the nine-label strategy and an independent POVM check."""

from functools import lru_cache
from itertools import product
import math

import numpy as np


QUESTIONS_PER_GAME = 243
CLASSICAL_LOSS = 1 / 27
TRACIAL_LOSS_FLOOR = (3 - 2 * math.sqrt(2)) / QUESTIONS_PER_GAME
FLAG = "flag"
PLUS_TRIPLE = (1, 1, 1)
EQUAL_SIGNS = ((1, 1), (-1, 1), (1, -1), (-1, -1))


def block_coordinates(rank, copy):
    """Indices of f, ++, -+, +-, -- for each block in one copy."""
    block = np.arange(2 * rank)
    flag = 9 * block + copy
    offsets = (copy + 1 + 2 * np.arange(4) + (block % 2)[:, None]) % 9
    data = 9 * (block // 2)[:, None] + offsets
    return np.column_stack((flag, data))


def strategy_loss(probabilities):
    """Exact loss formula for the explicit strategy, without shared ancillas."""
    rank = (len(probabilities) + 8) // 9
    coefficients = np.pad(np.sqrt(probabilities), (0, 18 * rank - len(probabilities)))
    total = 0.0
    for copy in range(9):
        flag, pp, mp, pm, mm = coefficients[block_coordinates(rank, copy)].T
        mean = (pp + mp + pm + mm) / 4
        q1 = (pp - mp + pm - mm) / 2
        q2 = (pp + mp - pm - mm) / 2
        q3 = (pp - mp - pm + mm) / 2
        total += np.sum((math.sqrt(2) * flag - mean) ** 2
                        + 8 * (q1 ** 2 + q3 ** 2) + (25 / 4) * q2 ** 2)
    return float(total / QUESTIONS_PER_GAME)


def triples(parity):
    return [signs for signs in product((1, -1), repeat=3) if math.prod(signs) == parity]


def projector(vector):
    return np.outer(vector, vector)


def spectral_projectors(observables, parity):
    outcomes = {}
    for signs in triples(parity):
        effect = np.eye(4)
        for sign, observable in zip(signs, observables):
            effect = effect @ ((np.eye(4) + sign * observable) / 2)
        outcomes[signs] = np.pad(effect, ((1, 0), (1, 0)))
    outcomes[FLAG] = np.diag([1., 0., 0., 0., 0.])
    return outcomes


@lru_cache(maxsize=1)
def measurements():
    """Construct the flagged magic-square projectors from the Pauli square."""
    identity = np.eye(2)
    x = np.array([[0., 1.], [1., 0.]])
    y = np.array([[0., -1j], [1j, 0.]])
    z = np.diag([1., -1.])
    square = [
        [np.kron(z, identity), np.kron(identity, z), np.kron(z, z)],
        [np.kron(identity, x), np.kron(x, identity), np.kron(x, x)],
        [np.kron(z, x), np.kron(x, z), np.kron(y, y).real],
    ]
    alice = {("R", i): spectral_projectors(row, 1) for i, row in enumerate(square)}
    bob = {("C", j): spectral_projectors(col, parity)
           for j, (col, parity) in enumerate(zip(zip(*square), (1, 1, -1)))}

    flag = np.eye(5)[0]
    ket_zero, ket_one = np.eye(2)
    x_vectors = {sign: np.array([1., sign]) / math.sqrt(2) for sign in (1, -1)}
    for label, triple in enumerate((PLUS_TRIPLE, (1, -1, -1))):
        c = np.r_[0., np.kron(ket_zero, x_vectors[(-1) ** label])]
        r = np.r_[0., np.kron(ket_zero, (ket_zero, ket_one)[label])]
        for side, kind, matched, vector in (
            (alice, "T", bob[("C", 0)], c),
            (bob, "U", alice[("R", 0)], r),
        ):
            effects = {answer: effect for answer, effect in matched.items()
                       if answer not in (FLAG, triple)}
            effects.update({sign: projector((flag + sign * vector) / math.sqrt(2))
                            for sign in (1, -1)})
            side[(kind, label)] = effects

    basis = np.column_stack([flag] + [np.r_[0., np.kron(x_vectors[s], x_vectors[t])]
                                      for s, t in EQUAL_SIGNS])
    for side in (alice, bob):
        for question, effects in side.items():
            side[question] = {answer: basis.T @ effect @ basis
                              for answer, effect in effects.items()}
    return alice, bob


def wins(x, y, a, b):
    kind_x, i = x
    kind_y, j = y
    if (kind_x, kind_y) == ("R", "C"):
        if a == FLAG or b == FLAG:
            return a == b
        return a[j] == b[i]
    if (kind_x, kind_y, j) == ("T", "C", 0):
        return b in (FLAG, (1, (-1) ** i, (-1) ** i)) if isinstance(a, int) else a == b
    if (kind_x, kind_y, i) == ("R", "U", 0):
        return a in (FLAG, (1, (-1) ** j, (-1) ** j)) if isinstance(b, int) else a == b
    if (kind_x, kind_y) == ("T", "U") and isinstance(a, int) and isinstance(b, int):
        return a * b == (-1) ** (i * j)
    return True


def fms_loss_kernel():
    """Sum the 25 flagged magic-square loss operators in the diagonal state sector."""
    alice, bob = measurements()
    kernel = np.zeros((5, 5))
    for x, a_effects in alice.items():
        for y, b_effects in bob.items():
            for a, a_effect in a_effects.items():
                for b, b_effect in b_effects.items():
                    if not wins(x, y, a, b):
                        kernel += a_effect * b_effect
    return kernel


def classical_strategy_lost_pairs():
    """Directly score the deterministic strategy attaining the classical value."""
    alice, bob = measurements()
    loses = []
    for copy in range(9):
        if copy == 0:
            a = {x: FLAG if x[0] == "R" else 1 for x in alice}
            b = {y: FLAG if y[0] == "C" else 1 for y in bob}
        else:
            t = 1 if 0 in {(copy + 1) % 9, (copy + 2) % 9} else -1
            a = {("R", 0): (-1, 1, -1), ("R", 1): (1, t, t), ("R", 2): (1, t, t),
                 ("T", 0): (-1, 1, -1), ("T", 1): (-1, 1, -1)}
            b = {("C", 0): (-1, 1, -1), ("C", 1): (1, t, t), ("C", 2): (-1, t, t),
                 ("U", 0): (-1, 1, -1), ("U", 1): (-1, 1, -1)}
        for x, a_answer in a.items():
            for y, b_answer in b.items():
                assert a_answer in alice[x] and b_answer in bob[y]
                if not wins(x, y, a_answer, b_answer):
                    loses.append((copy, x, y))
        if (a[("R", 0)] == FLAG) != (copy == 0):
            loses.append((copy, ("R", 0), "Lab"))
        if (a[("R", 1)] == PLUS_TRIPLE) != (0 in {(copy + 1) % 9, (copy + 2) % 9}):
            loses.append((copy, ("R", 1), "Lab"))
    return loses


def direct_strategy_loss(probabilities):
    """Sum all 243 question-pair losses by constructing effects from the rule table.

    The small dense matrices here are intended for checks on short spectra.
    Each copy's complemented effects sum to the identity on the padded space.
    """
    rank = (len(probabilities) + 8) // 9
    size = 18 * rank
    coefficients = np.pad(np.sqrt(probabilities), (0, size - len(probabilities)))
    pair_weights = np.outer(coefficients, coefficients)
    alice, bob = measurements()
    labels = np.diag(np.arange(size) % 9)
    total = 0.0

    def extend(side, indices):
        outside = np.setdiff1d(np.arange(size), indices)
        extended = {}
        for question, effects in side.items():
            extended[question] = {}
            for index, (answer, effect) in enumerate(effects.items()):
                matrix = np.zeros((size, size))
                matrix[indices[:, :, None], indices[:, None, :]] = effect
                if index == 0:
                    matrix[outside, outside] = 1
                extended[question][answer] = matrix
        return extended

    for copy in range(9):
        indices = block_coordinates(rank, copy)
        a_questions = extend(alice, indices)
        b_questions = extend(bob, indices)
        for x, a_effects in a_questions.items():
            for y, b_effects in b_questions.items():
                for a, a_effect in a_effects.items():
                    for b, b_effect in b_effects.items():
                        if not wins(x, y, a, b):
                            total += np.sum(pair_weights * a_effect * b_effect)
        for row, wanted, good in ((0, FLAG, {copy}),
                                  (1, PLUS_TRIPLE, {(copy + 1) % 9, (copy + 2) % 9})):
            for answer, effect in a_questions[("R", row)].items():
                for label in range(9):
                    if (answer == wanted) != (label in good):
                        lab_effect = np.diag(np.diag(labels) == label)
                        total += np.sum(pair_weights * effect * lab_effect)
    return float(total / QUESTIONS_PER_GAME)
