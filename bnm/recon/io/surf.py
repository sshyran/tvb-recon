import glob
import nibabel
from ..algo.geom import vertex_normals


def write_brain_visa_surf(fname, v, f):
    vn = vertex_normals(v, f)
    with open(fname, 'w') as fd:
        fd.write('- %d\n' % len(vn))
        for (vx, vy, vz), (nx, ny, nz) in zip(v, vn):
            fd.write('%f %f %f %f %f %f\n' % (vx, vy, vz, nx, ny, nz))
        fd.write('- %d %d %d\n' % ((len(f),)*3))
        for i, j, k in f:
            fd.write('%d %d %d\n' % (i, j, k))


def convert_fs_to_brain_visa(fs_surf):
    v, f = nibabel.freesurfer.read_geometry(fs_surf)
    write_brain_visa_surf(fs_surf + '.tri', v, f)


def convert_bem_to_tri(self):
    surfs_glob = '%s/%s/bem/watershed/*_surface-low' % (self.subjects_dir, self.subject)
    for surf_name in glob.glob(surfs_glob):
        convert_fs_to_brain_visa(surf_name)