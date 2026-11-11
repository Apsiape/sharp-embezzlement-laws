# Numerics for the paper: embezzlement error E = H^2 of half-chain states, and related statistics.
import sys, math, json, numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-claude/6fbe6fdd-66e6-555c-8896-876210039d7a/scratchpad/agents/a1')
from lib import OpenChain, zetas_majorana, zetas_from_matvec, spectrum_from_modes, H_functional, chain_modes

def F1(t):  # antiderivative of ln t
    t = np.asarray(t, dtype=float); out = np.zeros_like(t); m = t > 0
    out[m] = t[m] * np.log(t[m]) - t[m]; return out
def F2(t):  # antiderivative of ln^2 t
    t = np.asarray(t, dtype=float); out = np.zeros_like(t); m = t > 0
    lt = np.log(t[m]); out[m] = t[m] * lt ** 2 - 2 * t[m] * lt + 2 * t[m]; return out

def stats(p):
    """p sorted decreasing, sum ~ 1. Returns entropy (nats), b=-ln p_max, E=H^2, Var(log2 T), Var(K)."""
    p = np.asarray(p, dtype=float); p = p[p > 0]
    j = np.arange(1, len(p) + 1, dtype=float)
    S = float(-np.sum(p * np.log(p)))
    lp = np.log(p); mK = float(np.sum(p * lp)); varK = float(np.sum(p * (lp - mK) ** 2))
    m1 = float(np.sum(p * (F1(j) - F1(j - 1)))); m2 = float(np.sum(p * (F2(j) - F2(j - 1))))
    varX = (m2 - m1 ** 2) / math.log(2) ** 2
    E = H_functional(p) / 2.0
    return dict(S=S, b=-math.log(p[0]), E=E, varX=varX, varK=varK, mass=float(p.sum()), n=int(len(p)))

def xx_open(N, k=90):
    return zetas_from_matvec(OpenChain(N).block_matvec(N // 2), N // 2, k=min(k, N // 2 - 2))
def ising_open(Lspins, k=60):
    N = 2 * Lspins
    return zetas_majorana(OpenChain(N).block_matvec(N // 2), N // 2, k=min(k, N // 4 - 2))
def xx_gapped(N, delta, k=45):
    return zetas_from_matvec(OpenChain(N, delta=delta).block_matvec(N // 2), N // 2, k=min(k, N // 2 - 2))
