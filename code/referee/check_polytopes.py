"""Referee checks for section E (tetrahedron 7 vs 1) and the cube counterexample to (U). Fresh code."""
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial import ConvexHull
rng = np.random.default_rng(0)

def slice_poly(V, E, p, n):
    pts = []
    for i,j in E:
        a, b = (V[i]-p)@n, (V[j]-p)@n
        if a*b < 0: pts.append(V[i] + a/(a-b)*(V[j]-V[i]))
    return np.array(pts)
def slice_centroid(V, E, p, n):
    P = slice_poly(V, E, p, n)
    if len(P) < 3: return None, 0
    e1 = np.cross(n, [1,0,0] if abs(n[0])<0.9 else [0,1,0]); e1/=np.linalg.norm(e1); e2=np.cross(n,e1)
    Q = np.c_[(P-p)@e1, (P-p)@e2]; o = np.argsort(np.arctan2(Q[:,1]-Q[:,1].mean(), Q[:,0]-Q[:,0].mean())); Q=Q[o]
    x,y=Q[:,0],Q[:,1]; x2,y2=np.roll(x,-1),np.roll(y,-1); cr=x*y2-x2*y; A=cr.sum()/2
    cx=((x+x2)*cr).sum()/(6*A); cy=((y+y2)*cr).sum()/(6*A)
    return np.array([cx,cy]), abs(A)
def bary_dirs(V, E, p, starts=4000):
    def F(ang):
        th,ph=ang; n=np.array([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)])
        c,A=slice_centroid(V,E,p,n); return c if c is not None else np.array([9.,9.])
    sols=[]
    for _ in range(starts):
        z=rng.normal(size=3); z/=np.linalg.norm(z)
        r=least_squares(F,[np.arccos(z[2]),np.arctan2(z[1],z[0])],xtol=1e-14,ftol=1e-14,gtol=1e-14)
        if np.linalg.norm(r.fun)<1e-10:
            th,ph=r.x; n=np.array([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)])
            if n[np.argmax(abs(n))]<0: n=-n
            if all(min(np.linalg.norm(n-s),np.linalg.norm(n+s))>1e-6 for s in sols): sols.append(n)
    return sols
E6=[(i,j) for i in range(4) for j in range(i+1,4)]
for name,V in [("regular",np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]],float)),
               ("doc's tetra",np.array([[0,0,0],[3,0.2,0],[0.7,2.5,0.1],[0.9,0.8,2.2]],float)),
               ("random",rng.normal(size=(4,3)))]:
    g=V.mean(0); nv=0.85*V[0]+0.05*(V[1]+V[2]+V[3])
    print(name, "| at centroid:", len(bary_dirs(V,E6,g)), "| near vertex 0:", len(bary_dirs(V,E6,nv)))

# cube [-1,1]^3, p on the z-axis: cut-off volume for planes through p (exact via clipping)
def cut_volume(p, n):
    C=np.array([[x,y,z] for x in (-1,1) for y in (-1,1) for z in (-1,1)],float)
    E=[(i,j) for i in range(8) for j in range(i+1,8) if np.sum(C[i]!=C[j])==1]
    pts=[c for c in C if (c-p)@n>=0]+list(slice_poly(C,E,p,n))
    return ConvexHull(np.array(pts)).volume if len(pts)>=4 else 0.0
p=np.array([0,0,0.4])
vals=[]
for _ in range(20000):
    n=rng.normal(size=3); n/=np.linalg.norm(n)
    try: vals.append((min(cut_volume(p,n),8-cut_volume(p,n)),n))
    except Exception: pass
vmin=min(v for v,_ in vals); near=[n for v,n in vals if v<vmin+1e-9]
print("cube, p=(0,0,0.4): min cut volume %.6f (4(1-h)=%.6f); #random normals attaining it: %d of %d" % (vmin,4*(1-0.4),len(near),len(vals)))
