import numpy as np, raster, planar_reduce as pr, json
sec = raster.make_section(raster.torus_sdf(), 3.3, n=200)
out = {}
def run(name, p0, tguess, sign, h=0.02, steps=300):
    r = pr.inplane_dirs(sec, np.array(p0, float), n=96)
    t0 = min(r, key=lambda t: min(abs(t-tguess), abs(t-tguess+np.pi), abs(t-tguess-np.pi)))
    d = pr.ufun(t0)
    P, T, why = pr.trace(sec, np.array(p0, float), t0 if sign > 0 else t0, h=h*sign, steps=steps,
                         stop=lambda p: abs(p[2]) > 1.0 or np.hypot(p[0], 0) > 3.0)
    print(name, "start", p0, "deg %.1f" % np.degrees(t0), "->", np.round(P[-1], 3), why, len(P))
    out[name] = dict(P=P.tolist(), T=T.tolist(), why=why)
# fold branch: near-horizontal direction just above the equator
run("fold_a", [1.0, 0, 0.02], np.radians(172), +1)
run("fold_b", [1.0, 0, 0.02], np.radians(172), -1)
run("mid_a", [1.0, 0, 0.02], np.radians(122), +1)
run("mid_b", [1.0, 0, 0.02], np.radians(122), -1)
# main branch from several points
for i, (x, z) in enumerate([(2.0, 0.9), (2.6, 0.6), (1.4, 0.6), (2.9, 0.2), (1.1, 0.3)]):
    run("main%d_a" % i, [x, 0, z], np.radians(60), +1)
    run("main%d_b" % i, [x, 0, z], np.radians(60), -1)
json.dump(out, open("torus_traces.json", "w"))
