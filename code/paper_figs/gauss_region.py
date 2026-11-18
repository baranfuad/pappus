"""Figure for Example 2.1: the region under the bell curve y = exp(-x^2), its centroid,
the distance R from the centroid to the axis, and the axis of rotation y = -1.

Run:  python3 gauss_region.py      ->  writes ../../figures/gauss_region.png
"""
import os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ybar = 1/(2*np.sqrt(2))                    # centroid height of the region
x = np.linspace(-3.2, 3.2, 800); y = np.exp(-x*x)
BLUE, GREEN, RED, GREY = "#5b5fd6", "#2e9e4f", "#d64545", "#777777"

fig, ax = plt.subplots(figsize=(7.5, 4.6))
ax.fill_between(x, 0, y, color=BLUE, alpha=0.18)
ax.plot(x, y, color=BLUE, lw=3)
ax.axhline(0, color=GREY, lw=1.2); ax.plot([0, 0], [-1.25, 1.2], color=GREY, lw=1.2)
ax.text(3.35, -0.02, "$x$", fontsize=15, va="center"); ax.text(0.0, 1.24, "$y$", fontsize=15, ha="center")
ax.text(1.0, 0.82, "$y=e^{-x^2}$", color=BLUE, fontsize=16)
ax.plot([0, 0], [-1, ybar], color=GREEN, lw=1.8, ls="--")
ax.plot(0, ybar, "o", color=GREEN, ms=11, zorder=5)
ax.annotate("centroid $(0,\\,\\frac{1}{2\\sqrt{2}})$", (0.05, ybar), (1.45, 0.48), color=GREEN, fontsize=14,
            arrowprops=dict(arrowstyle="-", color=GREEN, lw=1.5), va="center")
ax.text(0.12, -0.5, "$R = 1+\\frac{1}{2\\sqrt{2}}$", color=GREEN, fontsize=15, va="center")
ax.plot([-3.2, 3.2], [-1, -1], color=RED, lw=2.5, ls=(0, (5, 4)))
ax.text(1.55, -0.86, "axis: $y=-1$", color=RED, fontsize=15)
ax.set_xlim(-3.3, 3.6); ax.set_ylim(-1.15, 1.35); ax.set_axis_off()
fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figures", "gauss_region.png")
fig.savefig(out, dpi=250, bbox_inches="tight"); print("wrote", os.path.normpath(out))
