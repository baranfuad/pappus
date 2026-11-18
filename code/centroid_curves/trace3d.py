"""Trace a centroid curve of a convex body in 3D by following a branch of barycentric directions."""
import numpy as np
from scipy.optimize import least_squares
from smooth import frame

def solve_dir(section, p, u0):
    e1, e2 = frame(u0)
    def res(x):
        u = u0 + x[0]*e1 + x[1]*e2; u /= np.linalg.norm(u)
        A, c = section(p, u); r = c - p
        return [r@e1, r@e2]
    s = least_squares(res, [0, 0], xtol=1e-14, ftol=1e-14, gtol=1e-14)
    u = u0 + s.x[0]*e1 + s.x[1]*e2
    return u/np.linalg.norm(u), np.linalg.norm(s.fun)

def trace(section, p0, u0, h, steps, stop):
    p = np.array(p0, float); u, _ = solve_dir(section, p, np.array(u0, float))
    P = [p.copy()]; U = [u.copy()]
    for k in range(steps):
        q = p + 0.5*h*u
        um, e = solve_dir(section, q, u); um = um if um@u > 0 else -um
        q = p + h*um
        un, e = solve_dir(section, q, um); un = un if un@u > 0 else -un
        if e > 1e-7: return np.array(P), np.array(U), "solver failed"
        p, u = q, un; P.append(p.copy()); U.append(u.copy())
        if stop(p): return np.array(P), np.array(U), "stop"
    return np.array(P), np.array(U), "max steps"
