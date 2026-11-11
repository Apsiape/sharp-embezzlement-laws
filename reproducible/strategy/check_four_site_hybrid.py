"""Exact SymPy checks for the four-site XX tilted-game separation."""

import itertools
import json

import sympy as sp


def occupation_index(bits):
    return sum(bit << (len(bits) - i - 1) for i, bit in enumerate(bits))


def spin_hamiltonian():
    hamiltonian = sp.zeros(16)
    for bits in itertools.product((0, 1), repeat=4):
        source = occupation_index(bits)
        for site in range(3):
            if bits[site] != bits[site + 1]:
                destination = source ^ (3 << (2 - site))
                hamiltonian[destination, source] = -1
    return hamiltonian


def slater_ground_state():
    # The six amplitudes are in the order of pairs 12, 13, 14, 23, 24, 34.
    amplitudes = (sp.sqrt(5) / 10, sp.Rational(1, 2), sp.sqrt(5) / 5,
                  sp.sqrt(5) / 5, sp.Rational(1, 2), sp.sqrt(5) / 10)
    state = sp.zeros(16, 1)
    for pair, amplitude in zip(itertools.combinations(range(4), 2), amplitudes):
        state[occupation_index(tuple(int(site in pair) for site in range(4)))] = amplitude
    return state


def fock_basis(orbitals):
    # Columns are the zero-, one-, and two-particle orbital occupation states.
    return sp.Matrix([
        [1, 0, 0, 0],
        [0, orbitals[1, 1], orbitals[1, 0], 0],
        [0, orbitals[0, 1], orbitals[0, 0], 0],
        [0, 0, 0, orbitals.det()],
    ])


def schmidt_basis_changes():
    left_orbitals = sp.Matrix([[1, -1], [1, 1]]) / sp.sqrt(2)
    right_orbitals = sp.Matrix([[1, 1], [-1, 1]]) / sp.sqrt(2)
    left_bits, right_bits = sp.zeros(4), sp.zeros(4)
    for n1, n2 in itertools.product((0, 1), repeat=2):
        source = occupation_index((n1, n2))
        left_bits[occupation_index((1 - n1, n2)), source] = (-1) ** ((1 - n1) * n2)
        right_bits[occupation_index((n1, 1 - n2)), source] = 1
    return left_bits * fock_basis(left_orbitals).T, right_bits * fock_basis(right_orbitals).T


def simplify_matrix(matrix):
    return matrix.applyfunc(sp.simplify)


def main():
    hamiltonian = spin_hamiltonian()
    ground = slater_ground_state()
    x = sp.Symbol("x")
    expected_polynomial = x ** 4 * (x ** 2 - 1) * (x ** 2 - 5) * (x ** 4 - 3 * x ** 2 + 1) ** 2
    assert sp.expand(hamiltonian.charpoly(x).as_expr() - expected_polynomial) == 0
    assert (ground.T * ground)[0] == 1
    assert simplify_matrix(hamiltonian * ground + sp.sqrt(5) * ground) == sp.zeros(16, 1)

    left, right = schmidt_basis_changes()
    assert simplify_matrix(left.T * left) == sp.eye(4)
    assert simplify_matrix(right.T * right) == sp.eye(4)
    local_change = sp.kronecker_product(left, right)
    schmidt = simplify_matrix(local_change * ground)
    nu = sp.Rational(1, 2) + 1 / sp.sqrt(5)
    expected_schmidt = sp.zeros(16, 1)
    for index, amplitude in zip((0, 5, 10, 15), (nu, sp.sqrt(5) / 10, sp.sqrt(5) / 10, 1 - nu)):
        expected_schmidt[index] = amplitude
    assert schmidt in (expected_schmidt, -expected_schmidt)

    z = sp.diag(1, -1)
    x_pauli = sp.Matrix([[0, 1], [1, 0]])
    a0, a1 = (sp.kronecker_product(pauli, sp.eye(2)) for pauli in (z, x_pauli))
    b0, b1 = (sp.kronecker_product((sp.sqrt(5) * z + sign * x_pauli) / sp.sqrt(6), sp.eye(2))
              for sign in (1, -1))
    for observable in (a0, a1, b0, b1):
        assert simplify_matrix(observable ** 2) == sp.eye(4)

    rows = []
    root_six_lower = sp.Rational(2449489, 10 ** 6)
    root_five_upper = sp.Rational(2236068, 10 ** 6)
    assert root_six_lower ** 2 < 6 and root_five_upper ** 2 > 5
    for alpha, robustness in ((sp.Integer(1), sp.Rational(39, 10000)),
                              (sp.Rational(5, 6), sp.Rational(58, 10000))):
        bell = (alpha * sp.kronecker_product(a0, sp.eye(4))
                + sp.kronecker_product(a0, b0 + b1)
                + sp.kronecker_product(a1, b0 - b1))
        physical_bell = local_change.T * bell * local_change
        score = sp.simplify((ground.T * physical_bell * ground)[0])
        assert sp.simplify(score - (2 * sp.sqrt(6) + 2 * alpha) / sp.sqrt(5)) == 0
        assert (2 + alpha) ** 2 > 8
        ground_loss = 1 / sp.Integer(270) + sp.Rational(9, 10) * (4 + alpha - score) / (2 * (4 + alpha))
        flat_floor = sp.Rational(9, 10) / (4 + alpha)
        gap = sp.simplify(flat_floor - ground_loss)
        score_lower = (2 * root_six_lower + 2 * alpha) / root_five_upper
        gap_lower = sp.Rational(9, 10) * (score_lower - 2 - alpha) / (2 * (4 + alpha)) - sp.Rational(1, 270)
        assert gap_lower > robustness
        rows.append({
            "alpha": str(alpha),
            "score_exact": str(score),
            "ground_loss": float(ground_loss),
            "all_flat_floor": float(flat_floor),
            "loss_gap": float(gap),
            "certified_trace_distance_margin": float(robustness),
        })
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
