import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/6fbe6fdd-66e6-555c-8896-876210039d7a/scratchpad/paper")
import math, numpy as np
from numlib import *
rng = np.random.default_rng(7); L2 = math.log(2)
def Eform(p):                      # E = 2 - sqrt2 * sum_m sqrt(p_m p_ceil(m/2)),  also ||lam - D lam||^2
    lam = np.sqrt(p); m = np.arange(1, len(p) + 1); half = lam[(m + 1) // 2 - 1]
    return 2 - math.sqrt(2) * float(np.sum(lam * half)), float(np.sum((np.concatenate([lam, np.zeros(len(lam))]) - np.repeat(lam, 2) / math.sqrt(2)) ** 2))
def spectra():
    for d in (1, 2, 3, 7, 64, 1000): yield f"flat{d}", np.full(d, 1.0 / d)
    for d in (5, 100, 5000): w = 1.0 / np.arange(1, d + 1); yield f"harmonic{d}", w / w.sum()
    for q in (0.2, 0.5, 0.9, 0.99): w = q ** np.arange(0, 3000); w = w[w > 1e-16]; yield f"geom{q}", w / w.sum()
    for k in range(40):
        d = int(rng.integers(2, 400)); w = np.sort(rng.random(d) ** rng.uniform(0.5, 12))[::-1] + 1e-300; yield f"rand{k}", w / w.sum()
    for n in (1024, 131072): p = np.load(f"spec_ising_{n}.npy"); yield f"ising{n}", p / p.sum()
    z = np.load("zeta_ising_4096.npy"); p, _ = spectrum_from_modes(np.concatenate([z, z])); yield "xx8192", p / p.sum()
worst = dict(i=0, iii=0, win=0, iv=0, v=9, nine=0, c=0, lip=0, inv=0)
for name, p in spectra():
    s = stats(p); E = s["E"]; E2, E3 = Eform(p)
    assert abs(E - E2) < 1e-10 and abs(E - E3) < 1e-10, (name, E, E2, E3)
    assert 0 < E <= 2 - math.sqrt(2) + 1e-12
    worst["iii"] = max(worst["iii"], p[0] / ((6 + 4 * math.sqrt(2)) * E))
    cs = np.concatenate([[0], np.cumsum(p), np.full(len(p) + 2, p.sum())])
    worst["win"] = max(worst["win"], max((cs[2 * r] - cs[r]) for r in range(1, len(p) + 1)) / math.sqrt(E), p[0] / 2 / math.sqrt(E))
    worst["iv"] = max(worst["iv"], (L2 / (2 * (s["S"] + 1 + 3 * L2))) / math.sqrt(E))
    worst["v"] = min(worst["v"], E * (4 * s["varX"] + 1))
    for k in (2, 3, 9): worst["inv"] = max(worst["inv"], abs(stats(np.repeat(p, k) / k)["E"] - E))
    # block-of-nine averaging and the chain of inequalities in the proof of Theorem 1(c)
    q = np.concatenate([p, np.zeros((-len(p)) % 9)]); q9 = np.repeat(q.reshape(-1, 9).mean(axis=1), 9)
    dist2 = float(np.sum((np.sqrt(q) - np.sqrt(q9)) ** 2)); E9 = stats(q9)["E"]
    worst["nine"] = max(worst["nine"], dist2 / (9 * p[0]))
    worst["lip"] = max(worst["lip"], abs(math.sqrt(E9) - math.sqrt(E)) / (2 * math.sqrt(dist2)) if dist2 > 0 else 0)
    bound = (math.sqrt(2 / 243) * math.sqrt(E9) + math.sqrt(dist2)) ** 2          # upper bound on the loss on psi
    worst["c"] = max(worst["c"], bound / (149 * E))
    hq = q.reshape(-1, 9).sum(axis=1); assert abs(stats(hq)["E"] - E9) < 1e-10      # E(chi x Phi_9) = E(chi)
print({k: round(v, 6) for k, v in worst.items()})
print("all ratios 'i..c' must be <= 1, 'v' >= 1, 'inv' ~ 0")
print("constants: 6+4sqrt2 =", 6 + 4 * math.sqrt(2), " (sqrt(2/243)*(1+6*sqrt(a)) + 3*sqrt(a))^2 =", (math.sqrt(2 / 243) * (1 + 6 * math.sqrt(6 + 4 * math.sqrt(2))) + 3 * math.sqrt(6 + 4 * math.sqrt(2))) ** 2, " 1+3ln2 =", 1 + 3 * L2, " ln2/(2(1+3ln2)) =", L2 / (2 * (1 + 3 * L2)))
