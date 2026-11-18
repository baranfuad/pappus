"""Checks for Theorem B (ellipsoid). a=sqrt2, b=c=1... and a general triaxial one."""
import numpy as np, smooth
from scipy.integrate import solve_ivp

def resid(phi, p, u):
    A, c = smooth.section(phi, p, u); r = c - p; return np.linalg.norm(r - (r@u)*u), A

# 1) barycentric direction at random p equals M p (normal of level ellipsoid)
a, b, c = 2.0, 1.3, 0.8
phi = smooth.ellipsoid(a, b, c); M = np.diag([1/a**2, 1/b**2, 1/c**2])
rng = np.random.default_rng(1)
for _ in range(3):
    p = rng.uniform(-1, 1, 3)*np.array([a, b, c])*0.5
    u = M@p; u /= np.linalg.norm(u)
    print("p", np.round(p, 3), "residual |c-p| for u=Mp:", resid(phi, p, u)[0])

# 2) parabola y = k x^2 is a centroid curve of the ellipsoid with a^2 = 2 b^2 (c arbitrary)
a, b, c = np.sqrt(2), 1.0, 0.7
phi = smooth.ellipsoid(a, b, c); k = 0.6
# endpoints: where parabola meets boundary
from scipy.optimize import brentq
X = brentq(lambda x: x**2/a**2 + (k*x*x)**2/b**2 - 1, 0, a)
xs = np.linspace(-X, X, 2001)[1:-1]
worst = 0; As = []
for x in xs:
    p = np.array([x, k*x*x, 0]); t = np.array([1, 2*k*x, 0]); t /= np.linalg.norm(t)
    e, A = resid(phi, p, t); worst = max(worst, e); As.append(A)
As = np.array(As)
ds = np.sqrt(1 + (2*k*xs)**2)
integral = np.trapezoid(As*ds, xs)
print("parabola k=%.2f: max centroid residual %.2e" % (k, worst))
print("  ∫A ds = %.6f   vol = 4/3 π abc = %.6f" % (integral, 4/3*np.pi*a*b*c))

# 3) do the normal planes along the parabola cross inside the body? (sweep not injective)
# planes at x1,x2: intersection line; check whether it meets K.
def planes_meet_inside(x1, x2):
    p1 = np.array([x1, k*x1**2, 0]); n1 = np.array([1, 2*k*x1, 0])
    p2 = np.array([x2, k*x2**2, 0]); n2 = np.array([1, 2*k*x2, 0])
    # both normals in xy-plane -> intersection line is vertical (z direction) through point q in xy plane
    Mx = np.array([n1[:2], n2[:2]]); rhs = np.array([n1[:2]@p1[:2], n2[:2]@p2[:2]])
    q = np.linalg.solve(Mx, rhs)
    return q[0]**2/a**2 + q[1]**2/b**2 < 1, q
hit = [(x1, x2) for x1 in np.linspace(-X*0.99, X*0.99, 60) for x2 in np.linspace(-X*0.99, X*0.99, 60)
       if x1 < x2 - 1e-9 and planes_meet_inside(x1, x2)[0]]
print("  pairs of normal planes that cross inside K:", len(hit), "e.g.", hit[:2])
