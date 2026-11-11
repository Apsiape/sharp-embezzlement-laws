import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/6fbe6fdd-66e6-555c-8896-876210039d7a/scratchpad/paper")
import time, json, math, numpy as np
from numlib import *
mode = sys.argv[1]
out = []
if mode == 'floor':
    for d in [float(x) for x in sys.argv[2].split(',')]:
        xi = 1 / math.asinh(d / 2); prev = None
        N = 2 ** int(math.ceil(math.log2(max(256, 8 * xi))))
        while True:
            t0 = time.time(); z = xx_gapped(N, d); p, lm = spectrum_from_modes(z, rel_thr=1e-11, cap=12_000_000)
            s = stats(p); s.update(model='gapped', size=N, delta=d, xi=xi, t=time.time() - t0)
            print(f"delta={d:<7} xi={xi:9.2f} N={N:>7} S={s['S']:.5f} b={s['b']:.5f} E={s['E']:.6f} E*S={s['E']*s['S']:.4f} varX={s['varX']:.4f} t={s['t']:.1f}", flush=True)
            if prev is not None and abs(prev - s['E']) < 1e-6: out.append(s); break
            if N >= 2 ** 18: out.append(s); break
            prev = s['E']; N *= 2
        if d == 0.02: np.save('spec_gapped_0.02.npy', p[:200000])
    json.dump(out, open(f'data_gapped_floor_{sys.argv[2]}.json', 'w'))
else:  # E versus N at fixed delta
    d = float(sys.argv[2])
    for k in range(3, int(sys.argv[3]) + 1):
        N = 2 ** k; t0 = time.time()
        z = chain_modes(N, d) if N <= 256 else xx_gapped(N, d)
        p, lm = spectrum_from_modes(z, rel_thr=1e-11, cap=12_000_000); s = stats(p); s.update(model='gappedN', size=N, delta=d, t=time.time() - t0); out.append(s)
        print(f"delta={d} N={N:>7} S={s['S']:.5f} E={s['E']:.6f} t={s['t']:.1f}", flush=True)
    json.dump(out, open(f'data_gappedN_{d}.json', 'w'))
