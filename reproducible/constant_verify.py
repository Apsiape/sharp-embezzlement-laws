"""Exact arithmetic for the flagged-sector and gluing bounds."""

from fractions import Fraction
from math import ceil


def add(*rows):
    return tuple(map(sum, zip(*rows, strict=True)))


def scale(factor, row):
    return tuple(factor * entry for entry in row)


def dot(left, right):
    return sum(a * b for a, b in zip(left, right, strict=True))


dimension = 11
basis = tuple(
    tuple(int(index == coordinate) for coordinate in range(dimension))
    for index in range(dimension)
)
a, b, u13, c, d, _, _, _, _, d_a, d_b = basis
row_sum = (2,) * 9 + (0, 0)
flag_bound_a = add(row_sum, scale(2, u13), scale(2, d_a))
flag_bound_b = add(row_sum, scale(2, u13))

# The three nonzero word-rewriting cases, with the sector-restriction
# error already added to each.
alice_c = (
    add(scale(4, add(a, c)), scale(2, d_a)),
    add(flag_bound_a, scale(4, a), scale(2, d_a)),
    add(flag_bound_a, scale(8, a), scale(4, add(b, c)), scale(2, d_a)),
)
alice_d = (
    add(flag_bound_a, scale(4, c), scale(2, d_a)),
    add(scale(4, add(b, d)), scale(2, d_a)),
    add(flag_bound_a, scale(4, add(b, d)), scale(2, d_a)),
)
bob_b = (
    add(scale(4, add(a, b)), scale(2, d_b)),
    add(flag_bound_b, scale(4, a), scale(2, d_b)),
    add(flag_bound_b, scale(8, a), scale(4, add(b, c)), scale(2, d_b)),
)
bob_d = (
    add(flag_bound_b, scale(4, b), scale(2, d_b)),
    add(scale(4, add(c, d)), scale(2, d_b)),
    add(flag_bound_b, scale(4, add(c, d)), scale(2, d_b)),
)


def replace_flag_mismatches(row):
    return (row[0] + row[9] + row[10], *row[1:9])


word_square_norms = tuple(
    tuple(dot(row, row) for row in map(replace_flag_mismatches, comparison))
    for comparison in (alice_c, alice_d, bob_b, bob_d)
)
assert word_square_norms == (
    (52, 144, 304),
    (112, 36, 144),
    (52, 108, 252),
    (92, 36, 124),
)

delta_a = scale(2, a)
delta_b = scale(2, add(b, d_b))
delta_c = scale(2, add(c, d_a))
delta_d = add(scale(2, d), b, c, d_a, d_b)
weighted_rows = [
    (Fraction(1), delta_a),
    (Fraction(7), delta_b),
    (Fraction(6), delta_c),
    (Fraction(13, 2), delta_d),
]
for coefficient, rows in (
    (Fraction(7, 24), bob_b),
    (Fraction(3, 10), alice_c),
    (Fraction(13, 24), alice_d),
    (Fraction(13, 20), bob_d),
):
    weighted_rows.extend((coefficient, row) for row in rows)

certificate = (1000, 623, 452, 634, 517, 226, 226, 226, 226, 833, 553)
diagonal = (Fraction(561, 1000),) + (1,) * 8 + (Fraction(439, 1000),) * 2
gram_certificate = add(
    *(scale(weight * dot(row, certificate), row) for weight, row in weighted_rows)
)
ratios = tuple(
    entry / (weight * coordinate)
    for entry, weight, coordinate in zip(
        gram_certificate, diagonal, certificate, strict=True
    )
)
assert all(ratio < Fraction("521.6") for ratio in ratios)
assert tuple(ceil(1000 * ratio) for ratio in ratios) == (
    521543, 521579, 521308, 521552, 521430, 521308,
    521308, 521308, 521308, 521546, 521180,
)
eta_bound = Fraction("11.42")
assert Fraction("521.6") < 4 * eta_bound**2

first_root = Fraction("11.464")
second_root = Fraction("11.594")
twice_sqrt_two = Fraction("2.828428")
assert first_root**2 > 1 + eta_bound**2
assert second_root**2 > 4 + eta_bound**2
assert twice_sqrt_two**2 > 8
magic_coefficient = first_root + twice_sqrt_two * (
    1 + 4 * first_root + 2 * second_root
)
assert magic_coefficient < Fraction("209.6")
assert Fraction("209.6")**2 + 160 < 210**2
assert Fraction(249, 4) + 35 * twice_sqrt_two / 2 < 112
global_constant = Fraction(243, 8) * (210**2 + 4 * 112)
assert global_constant == Fraction("1353145.5")

print("Word squared norms", word_square_norms)
print("Gram comparison ratios", tuple(round(float(ratio), 6) for ratio in ratios))
print("Magic coefficient upper bound", float(magic_coefficient))
print("Global constant", float(global_constant))
