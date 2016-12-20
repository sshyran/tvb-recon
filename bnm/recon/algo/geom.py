import numpy


# TODO test by generating points on unit sphere: vtx pos should equal normal
def vertex_normals(v, f):
    vf = v[f]
    fn = numpy.cross(vf[:, 1] - vf[:, 0], vf[:, 2] - vf[:, 0])
    vf = [set() for _ in range(len(v))]
    for i, fi in enumerate(f):
        for j in fi:
            vf[j].add(i)
    vn = numpy.zeros_like(v)
    for i, fi in enumerate(vf):
        fni = fn[list(fi)]
        norm = fni.sum(axis=0)
        norm /= numpy.sqrt((norm ** 2).sum())
        vn[i] = norm
    return vn


def merge_lh_rh(lv, lf, rv, rf, lrm, rrm):
    "Merge left and right hemisphere surfaces, and their region maps."
    v = numpy.r_[lv, rv]
    f = numpy.r_[lf, rf + lf.max()]
    rm = numpy.r_[lrm, rrm + lf.max()]
    return v, f, rm
