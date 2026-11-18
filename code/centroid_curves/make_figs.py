import numpy as np, json, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
# Figure 1: torus meridian plane
d = json.load(open("torus_exact_traces.json"))
fig, ax = plt.subplots(figsize=(8, 4.6))
th = np.linspace(0, 2*np.pi, 300)
for cx in (2, -2):
    ax.fill(cx + np.cos(th), np.sin(th), color="#dde6f2", zorder=0)
    ax.plot(cx + np.cos(th), np.sin(th), color="#5a6f8c", lw=1)
fx = [x for x, z in d["fold"] if z is not None]; fz = [z for x, z in d["fold"] if z is not None]
fx = np.r_[0, fx]; fz = np.r_[0, fz]
for s in (1, -1):
    ax.plot(fx, s*fz, color="#c0392b", lw=2, label="fold surface (3 directions inside, 1 outside)" if s == 1 else None)
    ax.plot(-fx, s*fz, color="#c0392b", lw=2)
ax.axhline(0, color="#888", lw=0.8, ls="--"); ax.axvline(0, color="#888", lw=0.8, ls="--")
for k, v in d["traces"].items():
    P = np.array(v["P"])
    col = "#1f6f4a" if k.startswith("main") else "#7d3c98"
    ax.plot(P[:, 0], P[:, 2], color=col, lw=1.3)
    if "fold" in v["why"]: ax.plot(P[-1, 0], P[-1, 2], "o", color="#c0392b", ms=4)
ax.set_xlim(-0.3, 3.2); ax.set_ylim(-1.1, 1.1); ax.set_aspect("equal")
ax.set_xlabel("distance from the torus axis (x)"); ax.set_ylabel("height z")
ax.set_title("Solid torus R=2, r=1, meridian plane: centroid curves end on a fold (red dots)")
ax.legend(loc="lower right", fontsize=8)
fig.tight_layout(); fig.savefig("figs/torus_meridian.png", dpi=160)
# zoom
ax.set_xlim(0.3, 1.9); ax.set_ylim(-0.25, 0.25); ax.set_aspect("auto"); ax.set_title("zoom on the fold region (vertical scale stretched)")
fig.savefig("figs/torus_fold_zoom.png", dpi=160)
# Figure 2: ellipsoid a=sqrt2, b=1 slice z=0: parabola centroid curves and the level ellipses
a, b = np.sqrt(2), 1
fig, ax = plt.subplots(figsize=(6, 4.4))
ax.plot(a*np.cos(th), b*np.sin(th), color="#5a6f8c")
for s in (0.25, 0.5, 0.75): ax.plot(s*a*np.cos(th), s*b*np.sin(th), color="#bbb", lw=0.7)
for k in (-1.2, -0.6, -0.25, 0, 0.25, 0.6, 1.2):
    x = np.linspace(-a, a, 2000); y = k*x*x; m = x**2/a**2 + y**2/b**2 <= 1
    ax.plot(x[m], y[m], color="#1f6f4a")
ax.plot([0, 0], [-b, b], color="#1f6f4a")
ax.set_aspect("equal"); ax.set_title("Ellipsoid with a = √2·b: every parabola y = kx² is a centroid curve")
fig.tight_layout(); fig.savefig("figs/ellipsoid_parabolas.png", dpi=160)
print("ok")
