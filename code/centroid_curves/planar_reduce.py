"""For bodies symmetric under y -> -y: barycentric directions u=(cos t,0,sin t) lying in the plane y=0,
at points p=(x,0,z). They are zeros of f(t) = (c(p,u)-p).(-sin t,0,cos t), with f(t+pi) = -f(t)."""
import numpy as np
from scipy.optimize import brentq

def ufun(t): return np.array([np.cos(t), 0.0, np.sin(t)])
def wfun(t): return np.array([-np.sin(t), 0.0, np.cos(t)])

def f(section, p, t):
    A, c = section(p, ufun(t)); return (c - p) @ wfun(t)

def inplane_dirs(section, p, n=72):
    ts = np.linspace(0, np.pi, n+1) + 0.0123456
    vals = [f(section, p, t) for t in ts]
    roots = []
    for i in range(n):
        if vals[i] == 0: roots.append(ts[i])
        elif vals[i]*vals[i+1] < 0:
            roots.append(brentq(lambda t: f(section, p, t), ts[i], ts[i+1], xtol=1e-10))
    return roots

def trace(section, p0, t0, h=0.02, steps=400, stop=lambda p: False):
    """follow the branch of directions continuously from (p0,t0); returns points, directions, and the reason it stopped."""
    p = np.array(p0, float); t = t0; P = [p.copy()]; T = [t]
    d = ufun(t)
    for k in range(steps):
        q = p + h*d
        # find root near previous angle
        lo, hi = t - 0.25, t + 0.25
        fl, fh = f(section, q, lo), f(section, q, hi)
        if fl*fh > 0:
            return np.array(P), np.array(T), "branch lost near %s" % np.round(q, 3)
        tn = brentq(lambda s: f(section, q, s), lo, hi, xtol=1e-10)
        dn = ufun(tn)
        if dn @ d < 0: dn = -dn
        p, t, d = q, tn, dn
        P.append(p.copy()); T.append(t)
        if stop(p): return np.array(P), np.array(T), "stop condition"
    return np.array(P), np.array(T), "max steps"
