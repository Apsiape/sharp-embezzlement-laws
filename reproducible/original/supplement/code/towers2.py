import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/6fbe6fdd-66e6-555c-8896-876210039d7a/scratchpad/paper")
import json, math, numpy as np
from numlib import *
import cl_and_pairs_lib as CLL
# gapped chain: are the one-body entanglement energies equally spaced?
for d, N in ((0.1, 2048), (0.01, 16384)):
    z = np.sort(xx_gapped(N, d))[::-1]; e = np.log((1 - z) / z)
    print("gapped delta", d, "eps_l/eps_0:", np.round(e[:8] / e[0], 5), "eps_0 =", round(float(e[0]), 5), flush=True)
z = np.sort(np.load('zeta_ising_131072.npy'))[::-1]; e = np.log((1 - z) / z)
print("ising 131072: pi^2/eps_0 =", math.pi ** 2 / e[0], " ln L =", math.log(131072), " kappa =", math.exp(math.pi ** 2 / e[0]) / 131072)
z = np.sort(np.load('zeta_ising_1024.npy'))[::-1]; e = np.log((1 - z) / z)
print("ising 1024: kappa =", math.exp(math.pi ** 2 / e[0]) / 1024)
def tower(eps, copies, lmax=400):
    e = (2 * np.arange(lmax) + 1) * eps; e = e[e < 40]
    z = 1 / (1 + np.exp(e)); z = np.concatenate([z] * copies)
    p, lm = spectrum_from_modes(z, rel_thr=1e-11, cap=12_000_000); s = stats(p); s['eps'] = eps; return s
out = {}
out['1'] = [tower(float(x), 1) for x in np.exp(np.linspace(math.log(5), math.log(0.25), 70))]
out['2'] = [tower(float(x), 2) for x in np.exp(np.linspace(math.log(5), math.log(0.40), 60))]
json.dump(out, open('data_towers.json', 'w'))
print("towers done", out['1'][-1]['b'], out['2'][-1]['b'], flush=True)
bs = sorted(set(list(np.round(np.arange(0.02, 1.4001, 0.02), 4)) + [1.5, 1.6, 1.7, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 8, 12, 16, 24, 32, 48, 64]))
cl = []
for b in bs:
    s = CLL.cl_stats(float(b)); cl.append(s)
json.dump(cl, open('data_cl_fine.json', 'w'))
print("cl done", len(cl))
