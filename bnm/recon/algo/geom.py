import numpy as np

def vertex_normals(v, f):
    vf = v[f]
    fn = np.cross(vf[:,1] - vf[:, 0], vf[:, 2] - vf[:, 0])
    vf = [set() for _ in range(len(v))]
    for i, fi in enumerate(f):
        for j in fi:
            vf[j].add(i)
    vn = np.zeros_like(v)
    for i, fi in enumerate(vf):
        fni = fn[list(fi)]
        norm = fni.sum(axis=0)
        norm /= np.sqrt((norm**2).sum())
        vn[i] = norm
    return vn
