"""Exact Schur-complement certificate for the universal Bell strategy bound."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class Quadratic:
    """A rational number plus a rational multiple of sqrt(2)."""

    rational: Fraction = Fraction(0)
    radical: Fraction = Fraction(0)

    def __add__(self, other: Quadratic) -> Quadratic:
        return Quadratic(self.rational + other.rational,
                         self.radical + other.radical)

    def __sub__(self, other: Quadratic) -> Quadratic:
        return Quadratic(self.rational - other.rational,
                         self.radical - other.radical)

    def __mul__(self, other: Quadratic) -> Quadratic:
        a, b = self.rational, self.radical
        c, d = other.rational, other.radical
        return Quadratic(a * c + 2 * b * d, a * d + b * c)

    def __truediv__(self, other: Quadratic) -> Quadratic:
        c, d = other.rational, other.radical
        denominator = c * c - 2 * d * d
        return self * Quadratic(c / denominator, -d / denominator)

    def positive(self) -> bool:
        a, b = self.rational, self.radical
        if a >= 0 and b >= 0:
            return a > 0 or b > 0
        if a <= 0 and b <= 0:
            return False
        if a > 0:
            return a * a > 2 * b * b
        return 2 * b * b > a * a


PIVOT_FLOORS = (3504, 6418, 10169, 8726, 11381, 10778, 11296, 10649, 11831,
                7230, 3906, 2185, 3944, 2287, 3512, 972, 4004, 2448, 145)


def verify() -> None:
    size = 19
    zero = Quadratic()
    inv_sqrt_two = Quadratic(radical=Fraction(1, 2))
    transform = [[Quadratic(Fraction(i == j)) for j in range(size)]
                 for i in range(size)]
    for i in range(size):
        transform[i][i // 2] -= inv_sqrt_two

    variance = [[Fraction(i == j) for j in range(size)] for i in range(size)]
    for block in (range(9), range(9, 18)):
        for i in block:
            for j in block:
                variance[i][j] -= Fraction(1, 9)
    variance[-1][-1] = Fraction(20, 9)

    scale = Quadratic(Fraction(15, 2))
    hessian = []
    for i in range(size):
        row = []
        for j in range(size):
            gram = sum((transform[k][i] * transform[k][j]
                        for k in range(size)), start=zero)
            row.append(scale * gram - Quadratic(variance[i][j]))
        hessian.append(row)

    for k, floor in enumerate(PIVOT_FLOORS):
        pivot = hessian[k][k]
        lower = Quadratic(Fraction(floor, 1000))
        upper = Quadratic(Fraction(floor + 1, 1000))
        assert (pivot - lower).positive() and (upper - pivot).positive()
        for i in range(k + 1, size):
            for j in range(k + 1, size):
                hessian[i][j] -= hessian[i][k] * hessian[k][j] / pivot


if __name__ == "__main__":
    verify()
    print("All 19 positive pivot intervals verified in Q(sqrt(2)).")
