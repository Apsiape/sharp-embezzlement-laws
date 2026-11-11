# Agent 1 numerics library: fast one-body entanglement spectra of free-fermion chains.
import math, numpy as np
from scipy.fft import dst
from scipy.sparse.linalg import LinearOperator, eigsh
from scipy.linalg import matmul_toeplitz

_SP = '/tmp/claude-0/-home-claude/6fbe6fdd-66e6-555c-8896-876210039d7a/scratchpad'
_src = open(_SP + '/tools/dilation_defect.py').read().split('print("reference spectra")')[0]
_ns = {}; exec(_src, _ns)
H_functional, spectrum_from_modes, chain_modes = _ns['H_functional'], _ns['spectrum_from_modes'], _ns['chain_modes']

def stats(p):
    p = p[p > 0]
    S = float(-np.sum(p * np.log(p)))            # nats
    lp = np.log(p); m = np.sum(p * lp)
    var = float(np.sum(p * (lp - m) ** 2))
    return dict(S=S, b=-math.log(p[0]), var=var, H=H_functional(p), mass=float(p.sum()), n=len(p))

# ---------- open chain, NN hopping, staggered potential delta, M particles, via DST-I ----------
class OpenChain:
    def __init__(self, N, M=None, delta=0.0):
        self.N = N; M = N // 2 if M is None else M
        k = np.arange(1, N + 1); q = np.pi * k / (N + 1); eps = -2 * np.cos(q)
        self.A = np.zeros(N); self.B = np.zeros(N)   # occupied orbital u = A_k phi_k + B_k phi_{N+1-k}
        if delta == 0.0:
            occ = np.argsort(eps, kind='stable')[:M]
            self.A[occ] = 1.0; self.pairs = None
        else:
            assert M == N // 2 and N % 2 == 0
            kk = np.arange(N // 2)                   # index of k (0-based), partner N-1-kk
            e = eps[kk]; r = np.sqrt(e ** 2 + delta ** 2)
            # block [[e, -d], [-d, -e]], lower eigvec ∝ (d, e + r) ... check: (e - (-r)) x - d y = 0 -> y = (e+r) x / d
            x = np.full_like(e, delta); y = -(e + r)  # dense code uses (-1)^(j-1) for 1-based j
            nrm = np.sqrt(x ** 2 + y ** 2); self.pairs = (kk, x / nrm, y / nrm)
        self.norm = math.sqrt(2.0 / (N + 1))
    def _dst(self, v):  # coefficients on normalized sine modes
        return dst(v, type=1) * (self.norm / 2)
    def matvec_full(self, v):
        w = self._dst(v)
        if self.pairs is None:
            c = self.A * w
        else:
            kk, a, b = self.pairs; part = self.N - 1 - kk
            s = a * w[kk] + b * w[part]
            c = np.zeros_like(w); c[kk] = a * s; c[part] = b * s
        return self._dst(c)
    def block_matvec(self, L):
        N = self.N
        def mv(x):
            x = np.asarray(x).ravel(); v = np.zeros(N); v[:L] = x
            return self.matvec_full(v)[:L]
        return mv
    def dense_block(self, L):
        mv = self.block_matvec(L); return np.array([mv(e) for e in np.eye(L)])

def toeplitz_sine_block(L, kF=np.pi / 2):
    n = np.arange(L, dtype=float); col = np.empty(L); col[0] = kF / np.pi
    col[1:] = np.sin(kF * n[1:]) / (np.pi * n[1:])
    return lambda x: matmul_toeplitz((col, col), np.asarray(x).ravel())

def zetas_from_matvec(mv, L, k=90, tol=1e-14):
    """Return one-body eigenvalues zeta (only the non-trivial ones) via eigsh on C(1-C)."""
    op = LinearOperator((L, L), matvec=lambda x: (lambda y: y - mv(y))(mv(x)), dtype=float)
    k = min(k, L - 2)
    x = eigsh(op, k=k, which='LA', tol=tol, return_eigenvectors=False)
    x = np.clip(x, 0, 0.25)
    return 0.5 - 0.5 * np.sqrt(1 - 4 * x)       # zeta <= 1/2 representative (spectrum symmetric)

def zetas_majorana(mv, L, k=60, tol=1e-14):
    """c=1/2 (TFI / Majorana chain) spectrum: nu = sing. values of K=(2C-1)[A_odd, A_even] (1-based sites)."""
    odd = np.arange(0, L, 2); even = np.arange(1, L, 2); m = len(odd)
    def kkT(u):
        v = np.zeros(L); v[odd] = u; w = 2 * mv(v) - v
        y = np.zeros(L); y[even] = w[even]; z = 2 * mv(y) - y
        return u - z[odd]
    op = LinearOperator((m, m), matvec=kkT, dtype=float)
    x = eigsh(op, k=min(k, m - 2), which='LA', tol=tol, return_eigenvectors=False)
    nu2 = np.clip(1 - x, 0, 1); return 0.5 - 0.5 * np.sqrt(nu2)
