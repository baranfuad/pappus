"""Run the interval-arithmetic certificates for the torus fold. Output -> torus_certificates.txt"""
import sys, time
from multiprocessing import Pool
from torus_arb import Mb, P, deg, arb, sign
def cert_box(args):
    x, z, a, b, psi0, psi1 = args
    psi = arb(0) if psi0 == psi1 == 0 else deg(psi0, psi1)
    stack, done = [(a, b)], []
    while stack:
        lo, hi = stack.pop(); v = Mb(P(x, z), deg(lo, hi), psi); s = sign(v)
        if s != "?": done.append((lo, hi, s, str(v))); continue
        if hi - lo < 1e-5: return (args, False, (lo, hi, str(v)))
        m = (lo+hi)/2; stack += [(m, hi), (lo, m)]
    return (args, True, done)
def cover(x, z, ranges, step, psi=(0, 0), pool=None):
    jobs = []
    for a, b in ranges:
        k = a
        while k < b - 1e-12: jobs.append((x, z, k, min(k+step, b), psi[0], psi[1])); k += step
    res = pool.map(cert_box, jobs, chunksize=4)
    ok = all(r[1] for r in res)
    signs = sorted([d for r in res if r[1] for d in r[2]])
    # summarise sign pattern as maximal runs
    runs = []
    for lo, hi, s, v in signs:
        if runs and runs[-1][2] == s and abs(runs[-1][1]-lo) < 1e-9: runs[-1][1] = hi
        else: runs.append([lo, hi, s])
    return ok, runs, len(signs), [r for r in res if not r[1]]
if __name__ == "__main__":
    out = open("torus_certificates.txt", "a"); T = time.time()
    def log(*a):
        s = " ".join(str(x) for x in a); print(s); out.write(s+"\n"); out.flush()
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 3) as pool:
        for spec in sys.argv[2:]:
            kind, x, z = spec.split(":")[0], float(spec.split(":")[1]), float(spec.split(":")[2])
            if kind == "pts":      # point evaluations (IVT sign changes)
                ts = [float(s) for s in spec.split(":")[3].split(",")]
                vals = pool.starmap(lambda_eval, [(x, z, t) for t in ts]) if False else [Mb(P(x, z), deg(t)) for t in ts]
                log("POINTS p=(%g,0,%g):" % (x, z), "; ".join("t=%g deg: Mb in %s [%s]" % (t, v, sign(v)) for t, v in zip(ts, vals)))
            elif kind == "cover":  # certify Mb != 0 on angle ranges, exclude given window
                rng = [tuple(map(float, r.split("-"))) for r in spec.split(":")[3].split(",")]
                ok, runs, nb, bad = cover(x, z, rng, float(spec.split(":")[4]), pool=pool)
                log("COVER p=(%g,0,%g) ranges %s: certified=%s, %d boxes, sign runs %s %s" % (x, z, rng, ok, nb, [(round(a,4), round(b,4), s) for a, b, s in runs], bad[:2]))
            elif kind == "patch":  # 2D patch of normals: t-range x psi-range
                rng = [tuple(map(float, spec.split(":")[3].split("-")))]; ps = tuple(map(float, spec.split(":")[5].split("/")))
                ok, runs, nb, bad = cover(x, z, rng, float(spec.split(":")[4]), psi=ps, pool=pool)
                log("PATCH p=(%g,0,%g) t in %s, psi in [%g,%g] deg: certified=%s, %d boxes, sign runs %s %s" % (x, z, rng, ps[0], ps[1], ok, nb, [(round(a,4), round(b,4), s) for a, b, s in runs], bad[:2]))
            log("  (elapsed %.0f s)" % (time.time()-T))
