"""Count barycentric directions in a (non-regular) tetrahedron: at the centroid and near a vertex."""
import numpy as np, box
V = np.array([[0, 0, 0], [3, 0.2, 0], [0.7, 2.5, 0.1], [0.9, 0.8, 2.2]], float)
E = [(V[i], V[j]) for i in range(4) for j in range(i+1, 4)]
g = V.mean(0)
def dirs_at(p, n=24):
    import types
    sols = box.bary_dirs(E, p, n=n)   # box.bary_dirs accepts an edge list
    return sols
for name, p in [("centroid", g), ("near vertex 0", 0.85*V[0] + 0.05*(V[1]+V[2]+V[3])), ("between", 0.5*g + 0.5*(0.7*V[0]+0.1*(V[1]+V[2]+V[3])))]:
    s = dirs_at(p); print(name, np.round(p, 3), len(s))
