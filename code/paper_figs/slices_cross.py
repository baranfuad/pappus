"""Figure for the paper (Section 6): the slices of the sqrt(2)-egg taken perpendicular to the
parabola y = k x^2 cross each other inside the body, yet their areas still add up to the volume.

We draw the picture in the plane z = 0. Each cutting plane is vertical (it contains the z-direction),
so it shows up as a straight line segment: the part of the line through (x, kx^2) with direction
(-2kx, 1) that lies inside the ellipse x^2/a^2 + y^2/b^2 <= 1.

Run:  python3 slices_cross.py      ->  writes ../../figures/slices_cross.png
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

a, b, k = np.sqrt(2), 1.0, 0.6
th = np.linspace(0, 2*np.pi, 400)
fig, ax = plt.subplots(figsize=(6, 4.4))
ax.fill(a*np.cos(th), b*np.sin(th), color="#dde6f2", zorder=0)
ax.plot(a*np.cos(th), b*np.sin(th), color="#5a6f8c", lw=1.2)

# where the parabola leaves the ellipse
X = max(r.real for r in np.roots([k*k/b**2, 0, 1/a**2, 0, -1]) if abs(r.imag) < 1e-12 and r.real > 0)
for x0 in np.linspace(-0.92*X, 0.92*X, 23):
    p = np.array([x0, k*x0*x0]); d = np.array([-2*k*x0, 1.0]); d /= np.linalg.norm(d)
    # solve |p + s d| on ellipse: quadratic in s
    A = d[0]**2/a**2 + d[1]**2/b**2; B = 2*(p[0]*d[0]/a**2 + p[1]*d[1]/b**2); C = p[0]**2/a**2 + p[1]**2/b**2 - 1
    s1, s2 = sorted(np.roots([A, B, C]).real)
    seg = np.array([p + s1*d, p + s2*d])
    ax.plot(seg[:, 0], seg[:, 1], color="#c0392b", lw=0.9, alpha=0.85)
x = np.linspace(-X, X, 400)
ax.plot(x, k*x*x, color="#1f6f4a", lw=2.2, label="centroid curve  $y=0.6\\,x^2$")
ax.plot([], [], color="#c0392b", lw=0.9, label="slices (seen edge-on)")
ax.set_aspect("equal"); ax.set_xlim(-1.55, 1.55); ax.set_ylim(-1.05, 1.05)
ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=8, frameon=False)
fig.tight_layout()
import os
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figures", "slices_cross.png")
fig.savefig(out, dpi=200)
print("wrote", os.path.normpath(out))
