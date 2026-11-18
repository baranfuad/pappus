"""Referee re-derivation of Lemma B1, Thm B2/B3, Cor B5 and Theorem A's sign conventions.
Independent of the original scripts (no shared code)."""
import sympy as sp, mpmath as mp, numpy as np

# ---- B1: centre of a planar section, computed directly by completing the square in the plane ----
a,b,c,t = sp.symbols('a b c t', positive=True)
u = sp.Matrix(sp.symbols('u1 u2 u3', real=True))
M = sp.diag(1/a**2,1/b**2,1/c**2)
# two vectors spanning the plane u.x = 0 (generic u)
U = sp.Matrix.hstack(sp.Matrix([u[1],-u[0],0]), sp.Matrix([u[2],0,-u[0]]))
p0 = t*u/(u.T*u)[0]                       # some point of the plane u.x = t
Q = U.T*M*U; g = U.T*M*p0
alpha = -Q.inv()*g                        # centre of the conic q(p0+U alpha) <= 1
centre = p0 + U*alpha
claim = t*M.inv()*u/(u.T*M.inv()*u)[0]
print("B1 centre == t M^-1 u/(u^T M^-1 u):", sp.simplify(centre-claim) == sp.zeros(3,1))

# ---- B2: p is the centre of the slice with normal Mp (and plane through p) ----
p = sp.Matrix(sp.symbols('p1 p2 p3', real=True)); uu = M*p; tt = (uu.T*p)[0]
print("B2 converse:", sp.simplify(tt*M.inv()*uu/(uu.T*M.inv()*uu)[0] - p) == sp.zeros(3,1))

# ---- B3: (x0 tau^(1/a^2), ...) has tangent parallel to M gamma ----
tau,x0,y0,z0 = sp.symbols('tau x0 y0 z0', positive=True)
gam = sp.Matrix([x0*tau**(1/a**2), y0*tau**(1/b**2), z0*tau**(1/c**2)])
print("B3 tangent x M gamma == 0:", sp.simplify(sp.diff(gam,tau).cross(M*gam)) == sp.zeros(3,1))

# ---- exact cap volume and slice area for an ellipsoid (via the map S = diag(a,b,c)) ----
mp.mp.dps = 30
def cap(n, q, A=(1,1,1)):
    """volume of {x in E : (x-q).n >= 0}"""
    S = mp.matrix([[A[0],0,0],[0,A[1],0],[0,0,A[2]]]); Sn = mp.matrix([A[0]*n[0],A[1]*n[1],A[2]*n[2]])
    nn = mp.sqrt(sum(x**2 for x in Sn)); d = (n[0]*q[0]+n[1]*q[1]+n[2]*q[2])/nn
    d = max(-1, min(1, d))
    return A[0]*A[1]*A[2]*mp.pi*(1-d)**2*(2+d)/3
def area(n, q, A=(1,1,1)):
    nn = mp.sqrt((A[0]*n[0])**2+(A[1]*n[1])**2+(A[2]*n[2])**2); d = (n[0]*q[0]+n[1]*q[1]+n[2]*q[2])/nn
    return A[0]*A[1]*A[2]*mp.pi*max(0,1-d**2)/nn
def centre(n, q, A=(1,1,1)):
    tt = n[0]*q[0]+n[1]*q[1]+n[2]*q[2]; Mi = [A[0]**2,A[1]**2,A[2]**2]
    den = sum(Mi[i]*n[i]**2 for i in range(3)); return [tt*Mi[i]*n[i]/den for i in range(3)]

# ---- Theorem A sign check on a NON-centroid curve (a helix-like curve), exact cap volumes ----
AX = (mp.mpf(2), mp.mpf('1.3'), mp.mpf('0.8'))
def gamma(s):  # not arc length; we reparametrise numerically
    return [mp.mpf('0.3')*mp.cos(s)+0.2*s, mp.mpf('0.25')*mp.sin(s), 0.1*s-0.05]
def frame(s):
    g1 = mp.diff(lambda x: gamma(x)[0], s); g2 = mp.diff(lambda x: gamma(x)[1], s); g3 = mp.diff(lambda x: gamma(x)[2], s)
    sp_ = mp.sqrt(g1**2+g2**2+g3**2); return [g1/sp_, g2/sp_, g3/sp_], sp_
def lhs_rhs(s0, s1):
    T0,_ = frame(s0); T1,_ = frame(s1)
    lhs = cap(T0, gamma(s0), AX) - cap(T1, gamma(s1), AX)
    def integrand(s):
        T, speed = frame(s)
        dT = [mp.diff(lambda x: frame(x)[0][i], s)/speed for i in range(3)]   # T' w.r.t. arc length
        q = gamma(s); A = area(T, q, AX); cc = centre(T, q, AX)
        # 1 - u kn + v kg  ==  1 - (c-gamma).T'  (frame-independent form)
        return A*(1 - sum((cc[i]-q[i])*dT[i] for i in range(3)))*speed
    return lhs, mp.quad(integrand, [s0, (s0+s1)/2, s1])
mp.mp.dps = 20
l, r = lhs_rhs(mp.mpf(-1), mp.mpf(1.5))
print("Theorem A on a non-centroid curve: LHS %s  RHS %s  diff %.1e" % (mp.nstr(l,12), mp.nstr(r,12), float(l-r)))

# ---- B5: sqrt2 parabola, exact areas, high-precision integral ----
mp.mp.dps = 30
for (aa,bb,cc_,k) in [(mp.sqrt(2),1,mp.mpf('0.7'),mp.mpf('0.6')), (mp.sqrt(2)*mp.mpf('1.1'),mp.mpf('1.1'),mp.mpf('0.5'),mp.mpf(3))]:
    AXp = (aa,bb,cc_)
    X = mp.findroot(lambda x: x**2/aa**2 + (k*x*x)**2/bb**2 - 1, aa*0.9)
    def res(x):
        q = [x, k*x*x, 0]; T = [1, 2*k*x, 0]; nT = mp.sqrt(1+4*k*k*x*x); T = [v/nT for v in T]
        cc2 = centre(T, q, AXp); return max(abs(cc2[i]-q[i]) for i in range(3))
    worst = max(res(x) for x in mp.linspace(-X*0.999, X*0.999, 201))
    I = mp.quad(lambda x: area([1/mp.sqrt(1+4*k*k*x*x), 2*k*x/mp.sqrt(1+4*k*k*x*x), 0], [x,k*x*x,0], AXp)*mp.sqrt(1+4*k*k*x*x), [-X,0,X])
    print("B5 parabola a=%s b=%s c=%s k=%s: max |centroid-p| %.1e ; ∫A ds - vol = %.1e" % (mp.nstr(aa,5),mp.nstr(bb,5),mp.nstr(cc_,5),k, float(worst), float(I-4*mp.pi*aa*bb*cc_/3)))

# ---- crossing-plane count quoted in section A (393 of 1770) ----
aa,bb,k = np.sqrt(2),1.0,0.6
from scipy.optimize import brentq
X = brentq(lambda x: x**2/aa**2+(k*x*x)**2/bb**2-1, 0, aa)
xs = np.linspace(-X*0.99, X*0.99, 60); cnt=0; tot=0
for i in range(60):
    for j in range(i+1,60):
        tot+=1; n1=np.array([1,2*k*xs[i]]); n2=np.array([1,2*k*xs[j]])
        qq=np.linalg.solve(np.array([n1,n2]), np.array([n1@[xs[i],k*xs[i]**2], n2@[xs[j],k*xs[j]**2]]))
        cnt += qq[0]**2/aa**2+qq[1]**2/bb**2<1
print("crossing pairs:", cnt, "of", tot)
