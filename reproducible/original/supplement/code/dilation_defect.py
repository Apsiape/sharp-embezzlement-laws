# Evaluate the dilation-defect functional of Theorem 1.4 of "Finite games that read the type
# of a von Neumann algebra",   H(h) = 2 * int_0^inf ( sqrt(p(t)) - sqrt(2 p(2t)) )^2 dt ,
# on the half-chain entanglement spectrum of free-fermion chains (my own exploratory numerics).
import numpy as np, math, sys

def H_functional(p):
    """p: Schmidt probabilities sorted decreasingly (type I: p(t) = p_j on (j-1, j])."""
    J = len(p)
    q = np.concatenate([p, np.zeros(2 * J + 2)])
    j = np.arange(1, J + 1)
    a = np.sqrt(q[j - 1])
    t1 = (a - np.sqrt(2 * q[2 * j - 2])) ** 2      # t in (j-1, j-1/2]: p(2t) = p_{2j-1}
    t2 = (a - np.sqrt(2 * q[2 * j - 1])) ** 2      # t in (j-1/2, j]:   p(2t) = p_{2j}
    return float(np.sum(t1 + t2))

def spectrum_from_modes(zeta, rel_thr=1e-11, cap=6_000_000):
    """Many-body Schmidt spectrum of a fermionic Gaussian state from one-body eigenvalues zeta."""
    zeta = np.clip(zeta, 1e-300, 1 - 1e-16)
    big = np.maximum(zeta, 1 - zeta); small = np.minimum(zeta, 1 - zeta)
    lam_max = float(np.exp(np.sum(np.log(big))))
    r = np.sort(small / big)[::-1]
    r = r[r > rel_thr]
    prods = np.array([1.0])
    for rk in r:
        new = prods * rk
        prods = np.concatenate([prods, new[new > rel_thr]])
        if len(prods) > cap: raise RuntimeError("too many states; raise rel_thr")
    p = np.sort(prods)[::-1] * lam_max
    return p, lam_max

def chain_modes(N, delta=0.0):
    """Open chain, nearest-neighbour hopping, staggered potential delta, half filling; left half."""
    h = np.zeros((N, N))
    i = np.arange(N - 1)
    h[i, i + 1] = h[i + 1, i] = -1.0
    h[np.arange(N), np.arange(N)] = delta * (-1.0) ** np.arange(N)
    w, v = np.linalg.eigh(h)
    occ = v[:, : N // 2]
    C = occ @ occ.T
    return np.linalg.eigvalsh(C[: N // 2, : N // 2])

def report(label, p, lam_max):
    S = float(-np.sum(p[p > 0] * np.log2(p[p > 0])))
    return f"{label:<26} states={len(p):>8} mass={p.sum():.9f} S={S:6.3f} ebits  -ln(lam_max)={-math.log(lam_max):6.3f}  H={H_functional(p):.5f}"

print("reference spectra")
for d in (2, 16, 256, 4096):
    p = np.full(d, 1.0 / d); print(report(f"flat, d={d}", p, 1.0 / d))
for d in (16, 256, 4096, 65536, 2 ** 20):
    w = 1.0 / np.arange(1, d + 1); p = w / w.sum(); print(report(f"1/t truncated, d={d}", p, p[0]), f" H*ln(d)={H_functional(p)*math.log(d):.3f}")
for n in (1, 2, 4, 8, 16, 32, 64):
    # n pairs with Schmidt weights (2/3, 1/3): multiplicities binom(n,k), weight (2/3)^(n-k) (1/3)^k
    ks = np.arange(n + 1)
    logw = (n - ks) * math.log(2 / 3) + ks * math.log(1 / 3)
    mult = np.array([math.comb(n, int(k)) for k in ks], dtype=float)
    if mult.sum() <= 4e6:
        p = np.repeat(np.exp(logw), mult.astype(int)); p = np.sort(p)[::-1]
        print(report(f"{n} pairs (2/3,1/3)", p, p[0]), f" H*n={H_functional(p)*n:.3f}")

print("\ncritical chain (delta = 0), left half of an open chain of N sites")
rows = []
for N in (8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096):
    p, lm = spectrum_from_modes(chain_modes(N))
    Hv = H_functional(p); rows.append((N, Hv))
    print(report(f"N={N}", p, lm), f" H*ln(N)={Hv*math.log(N):.3f}")
x = np.array([math.log(n) for n, _ in rows[3:]]); y = np.array([1 / h for _, h in rows[3:]])
a, b = np.polyfit(x, y, 1)
print(f"fit over N>=64: 1/H = {a:.3f} ln N + {b:.3f}")

print("\ngapped chain (staggered potential delta), N = 1024")
for delta in (0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.0):
    p, lm = spectrum_from_modes(chain_modes(1024, delta))
    print(report(f"delta={delta}", p, lm))
