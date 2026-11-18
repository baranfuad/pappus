"""Plane sections of a (possibly non-convex) body given by a signed distance function sdf(X) (<0 inside).
Area and centroid are computed on a grid in the plane with an anti-aliased indicator, which makes them
smooth functions of (p,u) up to O(h^2)."""
import numpy as np
from smooth import frame

def make_section(sdf, half_width, n=300):
    g = np.linspace(-half_width, half_width, n); h = g[1]-g[0]
    U, W = np.meshgrid(g, g, indexing='ij'); U = U.ravel(); W = W.ravel()
    def section(p, u):
        e1, e2 = frame(u)
        # centre the grid at the projection of the origin onto the plane (bodies are centred near 0)
        o = p - (p@e1)*e1 - (p@e2)*e2
        X = o + np.outer(U, e1) + np.outer(W, e2)
        w = np.clip(0.5 - sdf(X)/h, 0, 1)
        A = w.sum()*h*h
        if A == 0: return 0.0, p.copy()
        c = o + (w@U)*h*h/A*e1 + (w@W)*h*h/A*e2
        return A, c
    return section

def torus_sdf(R=2.0, r=1.0):
    def f(X):
        rho = np.sqrt(X[:, 0]**2 + X[:, 1]**2)
        return np.sqrt((rho-R)**2 + X[:, 2]**2) - r
    return f
