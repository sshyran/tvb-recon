import numpy
from bnm.recon.model.surface import Surface


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


# Merge left and right hemisphere surfaces, and their region maps.
def merge_lh_rh(lh_surface, rh_surface, left_region_mapping, right_region_mapping):
    out_surface = Surface([], [], [], None)
    out_surface.vertices = numpy.r_[lh_surface.vertices, rh_surface.vertices]
    out_surface.triangles = numpy.r_[lh_surface.triangles, rh_surface.triangles + lh_surface.vertices.max()]
    out_region_mapping = numpy.r_[left_region_mapping, right_region_mapping + lh_surface.triangles.max()]
    return out_surface, out_region_mapping
