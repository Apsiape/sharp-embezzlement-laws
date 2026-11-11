import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/6fbe6fdd-66e6-555c-8896-876210039d7a/scratchpad/paper")
import sys, time, json, numpy as np
from numlib import *
which = sys.argv[1]; kmin, kmax = int(sys.argv[2]), int(sys.argv[3])
out = []
for k in range(kmin, kmax + 1):
    t0 = time.time(); n = 2 ** k
    if which == 'xx':
        z = chain_modes(n) if n <= 256 else xx_open(n)
    else:
        if n <= 64:
            N = 2 * n; C = np.array([OpenChain(N).block_matvec(N // 2)(e) for e in np.eye(N // 2)])
            nu = np.linalg.svd((2 * C - np.eye(N // 2))[0::2, 1::2], compute_uv=False); z = (1 - nu) / 2
        else:
            z = ising_open(n)
    p, lm = spectrum_from_modes(z, rel_thr=1e-11, cap=12_000_000)
    s = stats(p); s.update(model=which, size=n, t=time.time() - t0); out.append(s)
    print(f"{which} size={n:>7} S={s['S']:.5f} b={s['b']:.5f} E={s['E']:.6f} E*S={s['E']*s['S']:.4f} varX={s['varX']:.4f} E*(4varX+1)={s['E']*(4*s['varX']+1):.4f} varK={s['varK']:.4f} n={s['n']} mass={s['mass']:.9f} t={s['t']:.1f}", flush=True)
    if k in (10, 17) : np.save(f'spec_{which}_{n}.npy', p[:200000])
json.dump(out, open(f'data_{which}_{kmin}_{kmax}.json', 'w'))
