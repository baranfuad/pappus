import numpy as np, json, warnings, torus_exact as te
warnings.filterwarnings("ignore")
# 1) fold curve z_f(x): largest z with 3 in-plane directions
fold = []
for x in np.linspace(0.1, 1.6, 16):
    lo, hi = 0.0, 0.3
    if len(te.dirs(np.array([x, 0, 1e-4]), n=120)) < 3: fold.append((x, None)); continue
    for _ in range(18):
        m = 0.5*(lo+hi)
        if len(te.dirs(np.array([x, 0, m]), n=120)) >= 3: lo = m
        else: hi = m
    fold.append((x, lo)); print("fold x=%.2f z=%.4f" % (x, lo), flush=True)
# 2) traces
out = {"fold": fold, "traces": {}}
def run(name, p0, tdeg, sgn, h=0.02, steps=200):
    p0 = np.array(p0, float); ds = te.dirs(p0, n=120)
    t0 = min(ds, key=lambda t: abs(((t - np.radians(tdeg)) + np.pi/2) % np.pi - np.pi/2))
    P, T, why = te.trace(p0, t0, sgn*h, steps)
    out["traces"][name] = dict(P=P.tolist(), T=T.tolist(), why=why)
    print(name, np.round(p0, 2), "%.1f deg" % np.degrees(t0), "->", np.round(P[-1], 3), why, flush=True)
for i, (x, z) in enumerate([(2.0, 0.95), (2.5, 0.8), (1.5, 0.8), (1.15, 0.5), (2.9, 0.3), (2.6, 0.2)]):
    run("main%d_down" % i, (x, 0, z), 60, -1)
    run("main%d_up" % i, (x, 0, z), 60, +1)
json.dump(out, open("torus_exact_traces.json", "w"))
