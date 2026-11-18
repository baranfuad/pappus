"""Exact (quadrature) version of the reduced problem for the solid torus (sqrt(x^2+y^2)-R)^2+z^2<=r^2.
Plane through p=(x,0,z) with normal u=(cos t,0,sin t). In-plane coordinates a (along w=(-sin t,0,cos t)) and b (along y).
For fixed a the slice in b has length L(a) in closed form; f(t) = centroid offset along w."""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
R, r = 2.0, 1.0

def L(a, p, t):
    qx = p[0] - a*np.sin(t); qz = p[2] + a*np.cos(t)
    if abs(qz) >= r: return 0.0
    s = np.sqrt(r*r - qz*qz)
    hi = (R+s)**2 - qx*qx
    if hi <= 0: return 0.0
    lo = max(0.0, (R-s)**2 - qx*qx)
    return 2*(np.sqrt(hi) - np.sqrt(lo))

def moments(p, t):
    c = np.cos(t)
    if abs(c) > 1e-9:
        a1, a2 = sorted([(-r - p[2])/c, (r - p[2])/c])
    else:
        a1, a2 = -(R+r)*2, (R+r)*2
    a1 = max(a1, -10); a2 = min(a2, 10)
    kw = dict(limit=400, epsabs=1e-13, epsrel=1e-12)
    m0 = quad(lambda a: L(a, p, t), a1, a2, **kw)[0]
    m1 = quad(lambda a: a*L(a, p, t), a1, a2, **kw)[0]
    return m0, m1

def f(p, t):
    m0, m1 = moments(p, t); return m1/m0 if m0 > 0 else 0.0

def area(p, t): return moments(p, t)[0]

def dirs(p, n=180):
    ts = np.linspace(0, np.pi, n+1) + 0.00731
    v = [f(p, s) for s in ts]
    return [brentq(lambda s: f(p, s), ts[i], ts[i+1], xtol=1e-13) for i in range(n) if v[i]*v[i+1] < 0]

def u(t): return np.array([np.cos(t), 0, np.sin(t)])

def trace(p0, t0, h, steps, zmax=1.0, win=0.2):
    p = np.array(p0, float); t = t0; d = u(t)*np.sign(h); P = [p.copy()]; T = [t]
    for k in range(steps):
        # midpoint (RK2) step on the branch
        def root_near(q, tg):
            lo, hi = tg-win, tg+win
            fl, fh = f(q, lo), f(q, hi)
            if fl*fh > 0: return None
            return brentq(lambda s: f(q, s), lo, hi, xtol=1e-13)
        q = p + 0.5*abs(h)*d
        tm = root_near(q, t)
        if tm is None: return np.array(P), np.array(T), "branch ends (fold) near %s" % np.round(q, 4)
        dm = u(tm); dm = dm if dm@d > 0 else -dm
        q = p + abs(h)*dm
        tn = root_near(q, tm)
        if tn is None: return np.array(P), np.array(T), "branch ends (fold) near %s" % np.round(q, 4)
        dn = u(tn); dn = dn if dn@d > 0 else -dn
        p, t, d = q, tn, dn; P.append(p.copy()); T.append(t)
        if abs(p[2]) > zmax or (np.hypot(p[0], 0)-R)**2 + p[2]**2 > r*r and abs(p[0]) > R: return np.array(P), np.array(T), "left body (outer side)"
    return np.array(P), np.array(T), "max steps"
