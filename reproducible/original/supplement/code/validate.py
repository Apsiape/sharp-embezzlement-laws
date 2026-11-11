# Independent validation of the free-fermion machinery against exact diagonalisation of spin chains.
import sys, math, numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-claude/6fbe6fdd-66e6-555c-8896-876210039d7a/scratchpad/agents/a1')
from lib import OpenChain, zetas_majorana, zetas_from_matvec, spectrum_from_modes, H_functional, chain_modes

sx = np.array([[0, 1], [1, 0.]]); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]); I2 = np.eye(2)
def op(o, j, L):
    out = np.array([[1.0 + 0j]])
    for k in range(L): out = np.kron(out, o if k == j else I2)
    return out
def schmidt(psi, L):
    M = psi.reshape(2 ** (L // 2), 2 ** (L // 2))
    s = np.linalg.svd(M, compute_uv=False) ** 2
    return np.sort(s)[::-1]

for L in (4, 8):
    # critical transverse-field Ising chain, open, L spins
    H = sum(-op(sx, j, L) @ op(sx, j + 1, L) for j in range(L - 1)) + sum(-op(sz, j, L) for j in range(L))
    w, v = np.linalg.eigh(H); p_ed = schmidt(v[:, 0], L)
    z = zetas_majorana(OpenChain(2 * L).block_matvec(L), L, k=min(60, L // 2 - 0)) if L >= 8 else None
    if z is None:
        # dense fallback for tiny sizes: singular values of the odd-even block of (2C-1) on the first L sites
        C = np.array([OpenChain(2 * L).block_matvec(L)(e) for e in np.eye(L)])
        K = (2 * C - np.eye(L))[0::2, 1::2]; nu = np.linalg.svd(K, compute_uv=False); z = (1 - nu) / 2
    p_ff, _ = spectrum_from_modes(z, rel_thr=1e-16)
    n = min(len(p_ed), len(p_ff))
    print(f"Ising L={L}: gap E1-E0={w[1]-w[0]:.4f}  max|p_ED - p_modes| over top {n} = {np.abs(p_ed[:n]-p_ff[:n]).max():.2e}  top3 ED {p_ed[:3]}")
    # XX chain, open, L sites
    H = sum(-0.5 * (op(sx, j, L) @ op(sx, j + 1, L) + op(sy, j, L) @ op(sy, j + 1, L)).real for j in range(L - 1))
    w, v = np.linalg.eigh(H); p_ed = schmidt(v[:, 0], L)
    p_ff, _ = spectrum_from_modes(chain_modes(L), rel_thr=1e-16)
    n = min(len(p_ed), len(p_ff))
    print(f"XX    L={L}: gap E1-E0={w[1]-w[0]:.4f}  max|p_ED - p_modes| over top {n} = {np.abs(p_ed[:n]-p_ff[:n]).max():.2e}  top3 ED {p_ed[:3]}")

# fast (DST + Lanczos) versus dense diagonalisation at N = 2048
N = 2048
zd = np.sort(np.minimum(chain_modes(N), 1 - chain_modes(N)))[::-1][:40]
zf = np.sort(zetas_from_matvec(OpenChain(N).block_matvec(N // 2), N // 2, k=60))[::-1][:40]
print("XX N=2048 fast vs dense one-body spectrum, max abs diff over top 40 modes:", np.abs(zd - zf).max())
pd_, _ = spectrum_from_modes(chain_modes(N)); pf, _ = spectrum_from_modes(zetas_from_matvec(OpenChain(N).block_matvec(N // 2), N // 2, k=60))
print("H^2 dense, fast:", H_functional(pd_) / 2, H_functional(pf) / 2)
