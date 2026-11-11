import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/6fbe6fdd-66e6-555c-8896-876210039d7a/scratchpad/paper")
import math, json, numpy as np
from scipy.special import i0e, i1e
from scipy.integrate import quad

# ---------- Calabrese-Lefevre spectrum: n(lambda) = I0(2 sqrt(b ln(lmax/lambda))), lmax = e^{-b}
xs = np.linspace(0, 400, 800001); logI0 = np.log(i0e(xs)) + xs          # ln I0(x), increasing
def x_of_lnt(lnt): return np.interp(lnt, logI0, xs)                      # inverse of ln I0
def cl_stats(b):
    lmax = math.exp(-b)
    p = lambda lnt: lmax * np.exp(-x_of_lnt(lnt) ** 2 / (4 * b))         # p(t) for t >= 1 as a function of ln t
    L = np.linspace(0, 380, 3_000_001); t = np.exp(L)
    pt = p(L); p2t = p(L + math.log(2))
    A = 0.5 * lmax                                                        # t in (0,1/2]
    Lh = np.linspace(math.log(0.5), 0, 200001)                            # t in (1/2,1]: p(t)=lmax, p(2t) smooth
    A += np.trapezoid(np.sqrt(lmax * p(Lh + math.log(2))) * np.exp(Lh), Lh)
    A += np.trapezoid(np.sqrt(pt * p2t) * t, L)                           # t >= 1
    E = 2 - 2 * math.sqrt(2) * A
    w = pt * t                                                            # density in ln t for t >= 1
    mass = lmax + np.trapezoid(w, L)
    S = lmax * b + np.trapezoid(w * (-np.log(np.maximum(pt, 1e-300))), L)
    # moments of ln T: T uniform on (0,1) with weight lmax, plus smooth part
    m1 = lmax * (-1.0) + np.trapezoid(w * L, L); m2 = lmax * 2.0 + np.trapezoid(w * L ** 2, L)
    varX = (m2 - m1 ** 2) / math.log(2) ** 2
    return dict(b=b, E=float(E), S=float(S), varX=float(varX), mass=float(mass))
out = []
for b in [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.25, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]:
    s = cl_stats(b); s['lognormal'] = 2 * (1 - math.exp(-math.log(2) ** 2 / (16 * b))); out.append(s)
    print(f"CL b={b:6.2f} mass={s['mass']:.6f} S={s['S']:.4f} (2b={2*b:.2f}) E={s['E']:.6f} E*S={s['E']*s['S']:.4f} E*8b/(ln2)^2={s['E']*8*b/math.log(2)**2:.4f} lognormal E={s['lognormal']:.6f} E(4varX+1)={s['E']*(4*s['varX']+1):.4f}")
json.dump(out, open('data_cl.json', 'w'))

# ---------- n copies of a two-qubit state with Schmidt weights (lam, 1-lam)
def E_blocks(logv, c):
    """logv[k]: log of the k-th distinct eigenvalue (decreasing); c[k]: cumulative count up to block k."""
    K = len(logv); A = 0.0; c0 = [0.0] + list(c)
    for k in range(K):
        a, bnd = 2 * c0[k], 2 * c0[k + 1]
        for m in range(k, K):
            lo, hi = max(a, c0[m]), min(bnd, c0[m + 1])
            if hi > lo:
                lg = 0.5 * (logv[k] + logv[m]) + math.log(hi - lo)
                if lg > -700: A += 0.5 * math.exp(lg)
            if c0[m + 1] >= bnd: break
    return 2 - 2 * math.sqrt(2) * A
pairs = []
for lam in (2 / 3, 0.8, 0.9, 0.75):
    for n in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1000):
        logv = [(n - k) * math.log(lam) + k * math.log(1 - lam) for k in range(n + 1)]
        c = np.cumsum([float(math.comb(n, k)) for k in range(n + 1)])
        E = E_blocks(logv, c); S = n * (-(lam * math.log(lam) + (1 - lam) * math.log(1 - lam)))
        pairs.append(dict(lam=lam, n=n, E=E, S=S))
    print(f"lam={lam:.4f}: " + "  ".join(f"n={d['n']}:{d['E']:.4f}" for d in pairs if d['lam'] == lam))
json.dump(pairs, open('data_pairs.json', 'w'))
# sanity: block routine against explicit enumeration for n = 8, lam = 2/3
from numlib import H_functional
n = 8; lam = 2 / 3
p = np.sort(np.concatenate([np.full(math.comb(n, k), lam ** (n - k) * (1 - lam) ** k) for k in range(n + 1)]))[::-1]
print("check n=8: explicit", H_functional(p) / 2, " blocks", [d['E'] for d in pairs if d['lam'] == lam and d['n'] == 8][0])
