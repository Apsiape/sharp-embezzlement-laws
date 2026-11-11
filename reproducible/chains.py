"""Dense and matrix-free spectra for the open XX and critical Ising chains."""

from dataclasses import dataclass
import math

import numpy as np
from scipy.fft import dst
from scipy.linalg import eigh_tridiagonal
from scipy.sparse.linalg import LinearOperator, eigsh

from spectrum import mode_probabilities, spectral_row


@dataclass(frozen=True)
class ModeSpectrum:
    occupation: np.ndarray
    full_count: int

    @property
    def omitted_mass_bound(self):
        """Bound conditional on the computed modes being the leading exact modes."""
        omitted = self.full_count - len(self.occupation)
        if not omitted:
            return 0.
        return -math.expm1(omitted * math.log1p(-float(np.min(self.occupation))))


class OpenChain:
    """Occupied-orbital projection for hopping -1 and alternating onsite potential."""

    def __init__(self, size, delta=0.):
        self.size = size
        self.half = size // 2
        self.delta = delta
        wave_numbers = np.pi * np.arange(1, self.half + 1) / (size + 1)
        energies = -2 * np.cos(wave_numbers)
        ratio = -delta / (np.hypot(energies, delta) - energies)
        self.a = 1 / np.sqrt(1 + ratio ** 2)
        self.b = self.a * ratio

    def block_matvec(self, vector):
        padded = np.zeros(self.size)
        padded[:self.half] = np.asarray(vector).ravel()
        sine = dst(padded, type=1, norm="ortho")
        low, high = sine[:self.half], sine[:self.half - 1:-1]
        if self.delta == 0:
            high.fill(0)
        else:
            occupied = self.a * low + self.b * high
            low[:] = self.a * occupied
            high[:] = self.b * occupied
        return dst(sine, type=1, norm="ortho")[:self.half]

    def dense_block(self):
        return np.column_stack([self.block_matvec(vector) for vector in np.eye(self.half)])


def dense_modes(size, delta=0.):
    onsite = delta * (-1.) ** np.arange(size)
    _, orbitals = eigh_tridiagonal(onsite, -np.ones(size - 1), select="i",
                                   select_range=(0, size // 2 - 1))
    orbitals = orbitals[:size // 2]
    occupations = np.linalg.eigvalsh(orbitals @ orbitals.T)
    return np.clip(np.minimum(occupations, 1 - occupations), 0., 0.5)


def leading_eigenvalues(matvec, dimension, count):
    operator = LinearOperator((dimension, dimension), matvec=matvec, dtype=float)
    start = np.random.default_rng(20261111).standard_normal(dimension)
    return eigsh(operator, k=count, which="LA", tol=1e-14, v0=start,
                 return_eigenvectors=False)


def xx_modes(size, delta=0., count=90, dense_cutoff=256):
    dimension = size // 2
    if size <= dense_cutoff or count >= dimension - 1:
        return ModeSpectrum(dense_modes(size, delta), dimension)
    block = OpenChain(size, delta).block_matvec

    def cross_product(vector):
        value = block(vector)
        return value - block(value)

    values = np.clip(leading_eigenvalues(cross_product, dimension, count), 0., 0.25)
    return ModeSpectrum(2 * values / (1 + np.sqrt(1 - 4 * values)), dimension)


def ising_modes(size, count=60, dense_cutoff=64):
    block = OpenChain(2 * size)
    dimension = size // 2
    if size <= dense_cutoff or count >= dimension - 1:
        covariance = 2 * block.dense_block() - np.eye(size)
        singular = np.linalg.svd(covariance[::2, 1::2], compute_uv=False)
        modes = np.clip((1 - singular) / 2, 0., 0.5)
        return ModeSpectrum(modes, dimension)

    def complement(vector):
        odd = np.zeros(size)
        odd[::2] = vector
        even = np.zeros(size)
        even[1::2] = 2 * block.block_matvec(odd)[1::2]
        return vector - 2 * block.block_matvec(even)[::2]

    values = np.clip(leading_eigenvalues(complement, dimension, count), 0., 1.)
    return ModeSpectrum(values / (2 * (1 + np.sqrt(1 - values))), dimension)


def chain_row(model, size, delta=0., threshold=1e-11, count=None):
    if size < 4 or size % 2 or (model == "xx-doubled" and size % 4):
        raise ValueError("Use an even size of at least 4, divisible by 4 for xx-doubled.")
    if model == "ising":
        modes = ising_modes(size, count or 60)
    elif model == "xx-doubled":
        half = ising_modes(size // 2, count or 60)
        modes = ModeSpectrum(np.tile(half.occupation, 2), 2 * half.full_count)
    elif model in ("xx", "gapped"):
        modes = xx_modes(size, delta, count or (45 if model == "gapped" else 90))
    else:
        raise ValueError(f"Unknown chain model {model!r}")
    return row_from_modes(modes, model, size, threshold, delta)


def row_from_modes(modes, model, size, threshold, delta=0.):
    probabilities = mode_probabilities(modes.occupation, threshold)
    row = spectral_row(probabilities, modes.omitted_mass_bound)
    row.update(model=model, size=size, modes_retained=len(modes.occupation), modes_total=modes.full_count)
    if model == "gapped":
        row.update(delta=delta, xi=1 / math.asinh(delta / 2))
    return row
