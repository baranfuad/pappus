"""Find all barycentric directions at p (critical directions of the cut-volume function)."""
import numpy as np, smooth
from scipy.optimize import least_squares

def sph(t, f): return np.array([np.sin(t)*np.cos(f), np.sin(t)*np.sin(f), np.cos(t)])
REF = np.array([0.31, 0.52, 0.79]); REF /= np.linalg.norm(REF)

def bary_dirs(section, p, n=10, tol=1e-9):
    sols = []
    def res(x):
        u = sph(*x); A, c = section(p, u); r = c - p
        e1 = np.cross(u, REF); e1 /= np.linalg.norm(e1); e2 = np.cross(u, e1)
        return [r@e1, r@e2]
    for t in np.linspace(0.1, np.pi/2, n//2):
        for f in np.linspace(0, 2*np.pi, n, endpoint=False):
            s = least_squares(res, [t, f], xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=60)
            if np.linalg.norm(s.fun) < tol:
                u = sph(*s.x); u = u if u @ REF > 0 else -u
                if all(min(np.linalg.norm(u-v), np.linalg.norm(u+v)) > 1e-5 for v in sols):
                    sols.append(u)
    return sols
