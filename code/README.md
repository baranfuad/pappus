# Code for "From a Bent Cylinder to Centroid Curves"

Every number and picture in the paper comes from one of these scripts. Appendix A of the paper has the same map with more explanation.

## Setup

```
pip install numpy scipy sympy mpmath matplotlib python-flint
```

Run each script **from inside its own folder** (they import each other by file name), for example:

```
cd code/centroid_curves
python3 sympy_checks.py
```

## Folders

- `centroid_curves/`: the original research code.
- `referee/`: an independent check, written from scratch without importing anything from `centroid_curves/`.
- `paper_figs/`: the figures made for the paper (each writes into `../../figures/`). Run them all with `cd code/paper_figs; for f in *.py; do python3 $f; done`.

## The four ideas the code uses

1. **Compute a slice.** Area and centroid of the cut of a body by a plane. Convex bodies: walk outwards from the point in 400 directions and add up thin triangles (`smooth.py`). Boxes and tetrahedra: exact polygons (`box.py`). Torus: add up straight chords whose lengths have a formula (`torus_exact.py`, `referee/torus_float.py`).
2. **Find barycentric directions.** For which cutting directions is the slice balanced exactly at the point? Root-finding from many starting guesses (`dirs.py`, `box.py`).
3. **Trace a centroid curve.** Small step along the current direction, find the nearest balanced direction at the new point, repeat. If none exists nearby, the curve has hit a fold (`trace3d.py`, `torus_exact.py`).
4. **Interval arithmetic.** The Arb library returns an interval guaranteed to contain the exact value, so a sign read from it is a proof (`referee/torus_arb.py`).

## Script map

| Script | What it does | Paper | Time |
|---|---|---|---|
| `centroid_curves/sympy_checks.py` | Symbolic: ellipsoid slice-centre formula; the corner triangle is balanced; x² − y² is constant on corner curves | Lemma 5.1, Prop. 8.1, Ex. 8.2 | seconds |
| `centroid_curves/ellipsoid_check.py` | Plane with normal Mp is balanced; the parabola y = 0.6x² in the √2-egg is a centroid curve and ∫A ds = vol; 393 crossing slice pairs | Thm 5.2, Ex. 5.7, Sec. 6 | ~1 min |
| `centroid_curves/smooth.py` | Library: slices of convex bodies; ellipsoid and superball | (library) | |
| `centroid_curves/dirs.py` | Library: all barycentric directions at a point | (library) | |
| `centroid_curves/trace3d.py` | Library: trace a centroid curve in 3-D | (library) | |
| `centroid_curves/superball_run.py` | Traces 5 curves of x⁴+y⁴+z⁴ ≤ 1 to the centre → `superball_traces.json` | Sec. 7 | minutes |
| `centroid_curves/box.py` | Library: exact polygon slices of boxes/polyhedra | (library) | |
| `centroid_curves/tetra_check.py` | Counts balanced planes in a tetrahedron (7 at centroid, 1 near a vertex) | Sec. 8 | minutes |
| `centroid_curves/raster.py` | Library: pixel-grid slices of non-convex bodies (first exploration) | (library) | |
| `centroid_curves/planar_reduce.py` | Library: directions in the symmetry plane y = 0 | (library) | |
| `centroid_curves/torus_trace.py` | First rough torus tracing → `torus_traces.json` (superseded) | exploration | minutes |
| `centroid_curves/torus_exact.py` | Library: high-precision torus slices, directions, tracing that stops at a fold | Sec. 9 | |
| `centroid_curves/torus_run.py` | Fold height at 16 x-values; traces from 6 surface points → `torus_exact_traces.json` | Sec. 9 | ~10 min |
| `centroid_curves/torus_more.py` | Traces the 2 curves starting inside the lens (purple) | Sec. 9 | minutes |
| `centroid_curves/make_figs.py` | Draws the torus and ellipsoid figures into `centroid_curves/figs/` | Figs. 1, 5 | seconds |
| `referee/check_ellipsoid.py` | Independent ellipsoid checks; Theorem 6.1 on a non-centroid curve (21 digits); √2-egg to 30 digits; 393/1770 | Sec. 5, 6 | seconds |
| `referee/check_polytopes.py` | Three tetrahedra (7 and 1); the cube counterexample to (U) | Sec. 7, 8 | ~10 min |
| `referee/check_superball.py` | Central-slice areas 3.708 / 4.410 / 4.443 and limiting directions | Sec. 7 | seconds |
| `referee/torus_float.py` | Independent torus slices; directions 40.8°, 46.3°, 113.3° and their disappearance | Sec. 9 | ~1 min |
| `referee/torus_arb.py` | Library: the torus integral in interval arithmetic | (library) | |
| `referee/torus_certify.py` | Produces the fold certificate → `torus_certificates.txt` | Sec. 9 | ~1 h |
| `paper_figs/slices_cross.py` | Draws Figure 4 (overlapping slices) into `figures/` | Fig. 4 | seconds |
| `paper_figs/gauss_region.py` | Bell-curve region with its centroid, the distance R and the axis y = −1 (Figure 3) | Ex. 2.1, Fig. 3 | seconds |

## Reproduce the main numbers

```
cd code/centroid_curves
python3 sympy_checks.py
python3 ellipsoid_check.py
python3 torus_run.py          # ~10 min
python3 torus_more.py
python3 make_figs.py

cd ../referee
python3 check_ellipsoid.py
python3 torus_float.py
python3 torus_certify.py 4 "pts:0.905:-0.098:30,43.5,80,130"
python3 torus_certify.py 4 "cover:0.905:-0.102:0-105,120-180:0.01"   # ~1 hour on 4 cores
```

`torus_certify.py` appends to `torus_certificates.txt`; the shipped file already holds the results quoted in the paper. The `.json` files are saved outputs, so `make_figs.py` works without re-running the torus traces.
