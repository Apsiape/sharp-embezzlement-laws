import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/6fbe6fdd-66e6-555c-8896-876210039d7a/scratchpad/paper")
import time, json, numpy as np
from numlib import *
# The open XX chain on N = 2L sites is two decoupled open critical Ising (Majorana) chains of L spins,
# so its one-body entanglement modes are the Ising modes, each taken twice.
out = []
for k in range(2, 18):
    t0 = time.time(); n = 2 ** k
    if n <= 64:
        N = 2 * n; C = np.array([OpenChain(N).block_matvec(N // 2)(e) for e in np.eye(N // 2)])
        nu = np.linalg.svd((2 * C - np.eye(N // 2))[0::2, 1::2], compute_uv=False); z = (1 - nu) / 2
    else:
        z = ising_open(n)
    np.save(f'zeta_ising_{n}.npy', z)
    p, lm = spectrum_from_modes(np.concatenate([z, z]), rel_thr=1e-11, cap=12_000_000)
    s = stats(p); s.update(model='xx', size=2 * n, t=time.time() - t0); out.append(s)
    print(f"xx(sq) size={2*n:>7} S={s['S']:.5f} b={s['b']:.5f} E={s['E']:.6f} varX={s['varX']:.4f} n={s['n']} mass={s['mass']:.9f} t={s['t']:.1f}", flush=True)
    json.dump(out, open('data_xx_sq.json', 'w'))
