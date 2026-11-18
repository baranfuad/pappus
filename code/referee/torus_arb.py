"""Computer-assisted (interval arithmetic) check of the torus fold, section F.
Uses Arb (python-flint): acb.integral returns a ball that is GUARANTEED to contain the exact integral,
with all rounding and truncation errors accounted for, provided the integrand is coded so that it
reports non-analyticity (we do: every sqrt(max(0,.)) returns NaN on a ball where its argument may vanish).
Integrand: same horizontal-chord formula as torus_float.py (meridian planes: d = e_y, m = (-sin t,0,cos t)).
Quantity: Mb(p,t) = int b L(b) db = first moment of the slice along the in-plane direction m.
Barycentric meridian direction  <=>  Mb = 0 (the other moment vanishes by the y -> -y symmetry)."""
from flint import arb, acb, ctx
import math, sys, time
ctx.prec = 80
NAN = acb(arb("nan"))
def psqrt(x, analytic):
    if analytic:
        if x.real > 0: return x.sqrt(analytic=True)
        if x.real < 0: return acb(0)
        return NAN
    r = x.real
    if r > 0: return acb(r.sqrt())
    if r < 0: return acb(0)
    return acb(arb(0).union(arb(r.upper()).sqrt()))
def nvec(t, psi):
    return (t.cos()*psi.cos(), psi.sin(), t.sin()*psi.cos())
def make_integrand(p, t, psi, weight):
    n1, n2, n3 = nvec(t, psi)
    if psi == 0:                         # exact meridian plane: d = e_y works for every t
        d1, d2 = arb(0), arb(1)
    else:
        s = (n1*n1 + n2*n2).sqrt(); d1, d2 = -n2/s, n1/s
    m1, m2, m3 = n2*0 - n3*d2, n3*d1 - n1*0, n1*d2 - n2*d1     # m = n x d (d3 = 0)
    def f(b, analytic):
        Q1, Q2, Q3 = p[0] + b*m1, p[1] + b*m2, p[2] + b*m3
        h = psqrt(1 - Q3*Q3, analytic)
        a0 = Q1*d1 + Q2*d2; rho2 = Q1*Q1 + Q2*Q2 - a0*a0
        Rp = psqrt((2+h)**2 - rho2, analytic); Rm = psqrt((2-h)**2 - rho2, analytic)
        L = 2*(Rp - Rm)
        return L*b if weight == "b" else L
    return f, float(m3.mid())
def Mb(p, t, psi=arb(0), weight="b"):
    f, m3 = make_integrand(p, t, psi, weight)
    if abs(m3) > 0.2:
        lo, hi = sorted([(-1-float(p[2].mid()))/m3, (1-float(p[2].mid()))/m3])
        lo, hi = max(lo-0.05, -8), min(hi+0.05, 8)
    else: lo, hi = -8, 8
    return acb.integral(f, lo, hi, rel_tol=2.0**-30, eval_limit=10**7, depth_limit=10**6).real
def deg(lo, hi=None):
    lo = arb(lo)*arb.pi()/180
    return lo if hi is None else arb(lo).union(arb(hi)*arb.pi()/180)
def sign(x): return "+" if x > 0 else "-" if x < 0 else "?"
def certify_nonzero(p, t0, t1, psi=arb(0), min_w=1e-4):
    """Cover [t0,t1] (degrees) by boxes on which Mb has a certified sign; returns (ok, list of (a,b,sign))."""
    stack, out = [(t0, t1)], []
    while stack:
        a, b = stack.pop(); v = Mb(p, deg(a, b), psi)
        s = sign(v)
        if s != "?": out.append((a, b, s)); continue
        if b - a < min_w: return False, (a, b, v)
        c = (a+b)/2; stack += [(c, b), (a, c)]
    out.sort(); return True, out
def P(x, z, y=0): return (arb(x), arb(y), arb(z))
if __name__ == "__main__":
    T = time.time()
    # sanity: slice by meridian plane t=0 (normal x) through (2,0,0) = two discs? no: plane x=2 cuts the tube
    # in a figure with area computed independently by the float code
    print("sanity A(p=(2,0,0), t=0) =", Mb(P(2,0), deg(0), weight="A"))
