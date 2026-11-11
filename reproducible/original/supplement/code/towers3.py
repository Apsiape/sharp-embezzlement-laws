import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/6fbe6fdd-66e6-555c-8896-876210039d7a/scratchpad/paper")
import json, math, numpy as np
from numlib import *
def tower(eps, copies, lmax=400):
    e = (2 * np.arange(lmax) + 1) * eps; e = e[e < 40]
    z = 1 / (1 + np.exp(e)); z = np.concatenate([z] * copies)
    p, lm = spectrum_from_modes(z, rel_thr=1e-11, cap=14_000_000); s = stats(p); s['eps'] = eps; return s
out = {'1': [], '2': []}
for copies, grid in ((1, (0.23, 0.21, 0.19, 0.175, 0.16, 0.148, 0.137, 0.127)), (2, (0.37, 0.34, 0.31, 0.285, 0.265, 0.25))):
    for eps in grid:
        try: s = tower(eps, copies)
        except RuntimeError as err: print("stop", copies, eps, err, flush=True); break
        out[str(copies)].append(s); json.dump(out, open('data_towers_ext.json', 'w'))
        print(f"copies={copies} eps={eps} b={s['b']:.4f} E={s['E']:.6f} 8bE/ln2^2={8*s['b']*s['E']/math.log(2)**2:.4f} n={s['n']} mass={s['mass']:.9f}", flush=True)
