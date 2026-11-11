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
