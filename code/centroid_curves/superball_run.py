import numpy as np, smooth, trace3d, json
phi = smooth.superball(4)
sec = lambda p, u: smooth.section(phi, p, u, R=2.0, it=42)
out = {}
rng = np.random.default_rng(5)
for k in range(5):
    p0 = rng.normal(size=3); p0 = np.abs(p0); p0 /= np.sum(p0**4)**0.25; p0 *= 0.97
    u0 = -p0**3; u0 /= np.linalg.norm(u0)          # inward normal of level set as initial guess
    P, U, why = trace3d.trace(sec, p0, u0, 0.01, 2000, stop=lambda p: np.linalg.norm(p) < 0.02)
    d = P[-1]/np.linalg.norm(P[-1])
    print("start", np.round(p0, 3), "end", np.round(P[-1], 4), why, "arrival direction", np.round(np.abs(d), 3), "tangent", np.round(np.abs(U[-1]), 3), flush=True)
    out[k] = dict(P=P.tolist(), why=why)
json.dump(out, open("superball_traces.json", "w"))
