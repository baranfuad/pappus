"""Symbolic checks: (1) ellipsoid section centroid formula; (2) barycentric plane at a cube corner."""
import sympy as sp
# (1) Ellipsoid: section of x^T M x <= 1 by u.x = t has centre t M^{-1}u/(u^T M^{-1} u).
a, b, c, t = sp.symbols('a b c t', positive=True)
u1, u2, u3 = sp.symbols('u1 u2 u3', real=True)
Minv = sp.diag(a**2, b**2, c**2); u = sp.Matrix([u1, u2, u3])
x0 = t*Minv*u/(u.T*Minv*u)[0]
# centre of a planar section of a quadric = point of the plane where grad q is parallel to u:
lam = sp.symbols('lam')
grad = sp.diag(1/a**2, 1/b**2, 1/c**2)*x0
print("(1) M x0 parallel to u:", sp.simplify(grad.cross(u)) == sp.zeros(3, 1), "; x0 on plane:", sp.simplify((u.T*x0)[0] - t) == 0)
# (2) Corner {x,y,z>=0}: plane through p with normal (1/p1,1/p2,1/p3) cuts triangle with centroid p.
p1, p2, p3 = sp.symbols('p1 p2 p3', positive=True)
n = sp.Matrix([1/p1, 1/p2, 1/p3]); d = n.dot(sp.Matrix([p1, p2, p3]))
V = [sp.Matrix([d/n[0], 0, 0]), sp.Matrix([0, d/n[1], 0]), sp.Matrix([0, 0, d/n[2]])]
cen = (V[0]+V[1]+V[2])/3
print("(2) corner triangle centroid - p:", sp.simplify(cen - sp.Matrix([p1, p2, p3])).T)
# (3) Gradient lines of xyz: x' = 1/x etc => x^2 - y^2 const
s = sp.symbols('s'); X = sp.Function('X')(s); Y = sp.Function('Y')(s)
print("(3) d/ds (x^2-y^2) along x'=1/x, y'=1/y:", sp.simplify(sp.diff(X**2 - Y**2, s).subs({sp.Derivative(X, s): 1/X, sp.Derivative(Y, s): 1/Y})))
