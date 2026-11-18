"""Referee's own float implementation for the solid torus (sqrt(x^2+y^2)-2)^2+z^2<=1.
General (non-horizontal) cutting planes, sliced into HORIZONTAL chords (derived independently of torus_exact.py).
Plane through p with unit normal n. d = horizontal unit vector in the plane, m = n x d.
Point X = p + a d + b m. For fixed b the chord is horizontal at height zeta = p3 + b m3,
its distance from the z-axis is rho, and it meets the torus in |a+a0| in [R-, R+],
R+- = sqrt(((2 +- h)^2 - rho^2)_+), h = sqrt((1-zeta^2)_+).
Area A = int L db, first moments M_b = int b L db, M_a = -int a0 L db, with L = 2(R+ - R-)."""
import numpy as np
from scipy.integrate import quad
def frame(n):
    n = np.asarray(n, float); n = n/np.linalg.norm(n)
    d = np.array([-n[1], n[0], 0.0]); d /= np.linalg.norm(d)
    if abs(n[1]) < 1e-15: d = np.array([0.0, 1.0, 0.0])
    m = np.cross(n, d); return n, d, m
def pieces(p, n, b):
    n, d, m = frame(n); Q = p + b*m; zeta = Q[2]
    h = np.sqrt(max(0.0, 1-zeta**2)); a0 = Q[:2]@d[:2]; rho2 = Q[:2]@Q[:2] - a0**2
    Rp = np.sqrt(max(0.0, (2+h)**2-rho2)); Rm = np.sqrt(max(0.0, (2-h)**2-rho2))
    return 2*(Rp-Rm), a0
def moments(p, n):
    n, d, m = frame(n); m3 = m[2]
    lo, hi = sorted([(-1-p[2])/m3, (1-p[2])/m3]); lo, hi = max(lo, -8), min(hi, 8)
    kw = dict(limit=500, epsabs=1e-14, epsrel=1e-13)
    A  = quad(lambda b: pieces(p, n, b)[0], lo, hi, **kw)[0]
    Mb = quad(lambda b: b*pieces(p, n, b)[0], lo, hi, **kw)[0]
    Ma = quad(lambda b: -pieces(p, n, b)[1]*pieces(p, n, b)[0], lo, hi, **kw)[0]
    return A, Ma, Mb
def nvec(t, psi=0.0): return np.array([np.cos(t)*np.cos(psi), np.sin(psi), np.sin(t)*np.cos(psi)])
def meridian_dirs(p, N=720):
    ts = np.linspace(0.0005, np.pi-0.0005, N); F = [moments(p, nvec(t))[2] for t in ts]
    roots = []
    for i in range(N-1):
        if F[i]*F[i+1] < 0:
            from scipy.optimize import brentq
            roots.append(np.degrees(brentq(lambda t: moments(p, nvec(t))[2], ts[i], ts[i+1], xtol=1e-12)))
    return roots, ts, np.array(F)
if __name__ == "__main__":
    for z in (-0.098, -0.100, -0.102):
        print("p=(0.905,0,%.3f):" % z, np.round(meridian_dirs(np.array([0.905, 0, z]), 360)[0], 3))
    # sanity: total volume via Cavalieri-type check  (area of slice through centre with normal z: annulus pi(9-1))
    print("A(z-normal at origin) =", moments(np.array([0,0,0.]), [0,0,1+1e-9]) if False else "skipped (vertical normal not supported)")
    print("A(x-normal through (2,0,0)) =", moments(np.array([2.,0,0]), [1,0,0])[0], " expected pi*(9-... ) see report")
