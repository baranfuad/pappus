"""Referee check of section D numbers for the superball |x|^4+|y|^4+|z|^4<=1 (central slices only)."""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
def A0(n):
    n = np.asarray(n, float); n /= np.linalg.norm(n)
    e1 = np.cross(n, [1,0,0] if abs(n[0]) < 0.9 else [0,1,0]); e1 /= np.linalg.norm(e1); e2 = np.cross(n, e1)
    r2 = lambda f: (np.sum((np.cos(f)*e1+np.sin(f)*e2)**4))**(-0.5)   # r(phi)^2
    return quad(lambda f: r2(f)/2, 0, 2*np.pi, limit=200)[0]
for n in ([1,0,0],[1,1,0],[1,1,1]): print("A0", n, round(A0(n), 4))
for e in ([1,0.3,0.1],[1,0.7,0.4]):
    e = np.array(e)/np.linalg.norm(e); best = None
    for s in range(30):
        x0 = np.random.default_rng(s).normal(size=3)
        r = minimize(lambda v: -A0(v)*(v@e)/np.linalg.norm(v), x0, method="Nelder-Mead", options=dict(xatol=1e-9, fatol=1e-12, maxiter=4000))
        if best is None or r.fun < best.fun: best = r
    v = best.x/np.linalg.norm(best.x); v *= np.sign(v@e); print("argmax A0(n)(n.e) for e", np.round(e,3), "->", np.round(v, 3))
