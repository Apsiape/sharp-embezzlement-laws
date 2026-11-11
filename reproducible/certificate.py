"""Recompute the exact 19-by-19 positive-definiteness certificate."""

from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class Quadratic:
    """Exact numbers a + b sqrt(2)."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __add__(self, other):
        return Quadratic(self.a + other.a, self.b + other.b)

    def __sub__(self, other):
        return Quadratic(self.a - other.a, self.b - other.b)

    def __mul__(self, other):
        return Quadratic(self.a * other.a + 2 * self.b * other.b,
                         self.a * other.b + self.b * other.a)

    def __truediv__(self, other):
        norm = other.a ** 2 - 2 * other.b ** 2
        return self * Quadratic(other.a / norm, -other.b / norm)

    def positive(self):
        if self.a >= 0 and self.b >= 0:
            return self.a > 0 or self.b > 0
        if self.a <= 0 and self.b <= 0:
            return False
        if self.a > 0:
            return self.a ** 2 > 2 * self.b ** 2
        return 2 * self.b ** 2 > self.a ** 2


PUBLISHED_FLOORS = (3504, 6418, 10169, 8726, 11381, 10778, 11296, 10649, 11831, 7230,
                    3906, 2185, 3944, 2287, 3512, 972, 4004, 2448, 145)


def pivots():
    size = 19
    t = [[Quadratic(Fraction(row == col), Fraction(-int(col == row // 2), 2))
          for col in range(size)] for row in range(size)]

    def q(row, col):
        if row == col == size - 1:
            return Fraction(20, 9)
        return Fraction(row == col) - Fraction(1, 9) if row // 9 == col // 9 < 2 else Fraction(0)

    scale = Quadratic(Fraction(15, 2))
    matrix = [[scale * sum((t[k][row] * t[k][col] for k in range(size)), Quadratic())
               - Quadratic(q(row, col)) for col in range(size)] for row in range(size)]

    for k in range(size):
        pivot = matrix[k][k]
        yield pivot
        for i in range(k + 1, size):
            factor = matrix[i][k] / pivot
            for j in range(k + 1, size):
                matrix[i][j] = matrix[i][j] - factor * matrix[k][j]


def verify_certificate():
    for pivot, expected in zip(pivots(), PUBLISHED_FLOORS, strict=True):
        assert (pivot - Quadratic(Fraction(expected, 1000))).positive()
        assert (Quadratic(Fraction(expected + 1, 1000)) - pivot).positive()
    return PUBLISHED_FLOORS


if __name__ == "__main__":
    print(verify_certificate())
