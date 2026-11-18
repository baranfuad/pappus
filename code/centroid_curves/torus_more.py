import numpy as np, json, warnings, torus_exact as te
warnings.filterwarnings("ignore")
d = json.load(open("torus_exact_traces.json"))
for name, tdeg in [("branch161", 161.5), ("branch125", 125.5)]:
    for sgn in (+1, -1):
        p0 = np.array([1.0, 0, 0.05]); ds = te.dirs(p0, n=180)
        t0 = min(ds, key=lambda t: abs(np.degrees(t)-tdeg))
        P, T, why = te.trace(p0, t0, sgn*0.01, 400)
        d["traces"]["%s_%+d" % (name, sgn)] = dict(P=P.tolist(), T=T.tolist(), why=why)
        print(name, sgn, "%.2f" % np.degrees(t0), "->", np.round(P[-1], 3), why, len(P), flush=True)
json.dump(d, open("torus_exact_traces.json", "w"))
