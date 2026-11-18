"""Exact plane sections of a box [-h1,h1]x[-h2,h2]x[-h3,h3] and barycentric directions."""
import numpy as np
from itertools import product
from scipy.optimize import least_squares

def box_edges(h):
    V = np.array(list(product(*[(-x, x) for x in h])), float)
    E = []
    for i in range(8):
        for j in range(i+1, 8):
            if np.sum(V[i] != V[j]) == 1: E.append((V[i], V[j]))
    return E

def section(h, p, u):
    """area and centroid of box ∩ {x : u.(x-p)=0}"""
    E = box_edges(h) if not isinstance(h, list) else h
    d = u @ p
    pts = []
    for a, b in E:
        fa, fb = u @ a - d, u @ b - d
        if fa == 0: pts.append(a)
        if fa * fb < 0: pts.append(a + (b - a) * fa / (fa - fb))
    if len(pts) < 3: return 0.0, p.copy()
    P = np.array(pts)
    m = P.mean(0)
    e1 = np.cross(u, [1, 0, 0] if abs(u[0]) < .9 else [0, 1, 0]); e1 /= np.linalg.norm(e1)
    e2 = np.cross(u, e1)
    ang = np.arctan2((P - m) @ e2, (P - m) @ e1)
    P = P[np.argsort(ang)]
    A = 0; C = np.zeros(3)
    for i in range(len(P)):
        a, b = P[i] - m, P[(i+1) % len(P)] - m
        t = 0.5 * np.linalg.norm(np.cross(a, b))
        A += t; C += t * (m + (a + b) / 3)
    return A, C / A

def sph(t, f): return np.array([np.sin(t)*np.cos(f), np.sin(t)*np.sin(f), np.cos(t)])

def F(h, p, u):
    A, c = section(h, p, u)
    r = c - p
    return r - (r @ u) * u

def bary_dirs(h, p, n=24, tol=1e-10):
    """all barycentric normal directions u (mod ±) at p, by multistart root finding."""
    E = box_edges(h) if not isinstance(h, list) else h
    sols = []
    for t in np.linspace(0.05, np.pi/2, n//2):
        for f in np.linspace(0, 2*np.pi, n, endpoint=False):
            def res(x):
                u = sph(*x); e1 = np.cross(u, [0.3, 0.5, 0.81]); e1/=np.linalg.norm(e1); e2=np.cross(u,e1)
                r = F(E, p, u); return [r @ e1, r @ e2]
            s = least_squares(res, [t, f], xtol=1e-14, ftol=1e-14, gtol=1e-14)
            if np.linalg.norm(s.fun) < tol:
                u = sph(*s.x); u = u if (u[2] > 1e-12 or (abs(u[2]) <= 1e-12 and (u[1] > 0 or (abs(u[1])<1e-12 and u[0]>0)))) else -u
                if all(min(np.linalg.norm(u - v), np.linalg.norm(u + v)) > 1e-6 for v in sols): sols.append(u)
    return sols
