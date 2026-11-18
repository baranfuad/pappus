"""Plane sections of a convex body K = {x : phi(x) <= 1} through an interior point p.
Uses polar coordinates about p in the plane: r(theta) found by bisection (valid since K convex, p interior)."""
import numpy as np

NT = 400
TH = np.linspace(0, 2*np.pi, NT, endpoint=False)

def frame(u):
    a = np.array([1., 0, 0]) if abs(u[0]) < 0.9 else np.array([0, 1., 0])
    e1 = np.cross(u, a); e1 /= np.linalg.norm(e1); e2 = np.cross(u, e1)
    return e1, e2

def section(phi, p, u, R=10.0, it=60):
    """area A and centroid c of K ∩ {u.(x-p)=0}; p must be inside K."""
    e1, e2 = frame(u)
    D = np.outer(np.cos(TH), e1) + np.outer(np.sin(TH), e2)   # NT x 3
    lo = np.zeros(NT); hi = np.full(NT, R)
    for _ in range(it):
        mid = 0.5*(lo+hi)
        inside = phi(p + mid[:, None]*D) <= 1
        lo = np.where(inside, mid, lo); hi = np.where(inside, hi, mid)
    r = 0.5*(lo+hi)
    dth = 2*np.pi/NT
    A = np.sum(r**2/2)*dth
    m = (np.sum((r**3/3)[:, None]*D, 0))*dth
    return A, p + m/A

def ellipsoid(a, b, c):
    s = np.array([a, b, c], float)
    return lambda X: np.sum((np.atleast_2d(X)/s)**2, axis=-1)

def superball(q=4, s=(1, 1, 1)):
    s = np.array(s, float)
    return lambda X: np.sum(np.abs(np.atleast_2d(X)/s)**q, axis=-1)
