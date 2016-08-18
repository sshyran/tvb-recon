import os
import sys
import glob
import numpy as np
import nibabel
import scipy.io
import warnings


try:
    import gdist
except ImportError:
    warnings.warn('Geodesic distance module unavailable; please pip install gdist.')


SUBJECTS_DIR, SUBJECT = [os.environ[key] for key in 'SUBJECTS_DIR SUBJECT'.split()]


def vertex_normals(v, f):
    vf = v[f]
    fn = np.cross(vf[:,1] - vf[:, 0], vf[:, 2] - vf[:, 0])
    vf = [set() for _ in xrange(len(v))]
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
    v, f = freesurfer.read_geometry(fs_surf)
    write_brain_visa_surf(fs_surf + '.tri', v, f)


def read_annot(hemi, annot_name):
    annot_fname = '%s.%s.annot' % (hemi, annot_name)
    annot_path = os.path.join(SUBJECTS_DIR, SUBJECT, 'label', annot_fname)
    return freesurfer.read_annot(annot_path)


def annot_to_lut(hemi, annot_name, lut_path):
    _, ctab, names = read_annot(hemi, annot_name)
    with open(lut_path, 'w') as fd:
        for name, (r, g, b, a, id) in zip(names, ctab):
            fd.write('%d\t%s\t%d %d %d %d\n' % (id, name, r, g, b, a))


def annot_to_conn_conf(hemi, annot_name, conn_conf_path):
    _, _, names = read_annot(hemi, annot_name)
    with open(conn_conf_path, 'w') as fd:
        for id, name in enumerate(names):
            fd.write('%d\t%s\n' % (id, name))


def compute_gdist_mat(surf_name='pial', max_distance=40.0):
    max_distance = float(max_distance) # in case passed from sys.argv
    for h in 'rl':
        surf_path = '%s/%s/surf/%sh.%s' % (SUBJECTS_DIR, SUBJECT, h, surf_name)
        v, f = nibabel.freesurfer.read_geometry(surf_path)
        mat_path = '%s/%s/surf/%sh.%s.gdist.mat' % (SUBJECTS_DIR, SUBJECT, h, surf_name)
        mat = gdist.local_gdist_matrix(v, f.astype('<i4'), max_distance=40.0)
        scipy.io.savemat(mat_path, {'gdist': mat})


def convert_bem_to_tri():
    surfs_glob = '%s/%s/bem/watershed/*_surface-low' % (SUBJECTS_DIR, SUBJECT)
    for surf_name in glob.glob(surfs_glob):
        utils.convert_fs_to_brain_visa(surf_name)


def gen_head_model():
    surfs_glob = '%s/%s/bem/watershed/*_surface-low.tri' % (SUBJECTS_DIR, SUBJECT)
    surfs = glob.glob(surfs_glob)

    if len(surfs) == 0:
        raise Exception('tri surfaces not found!')

    hm_base = '%s/%s/bem/head_model' % (SUBJECTS_DIR, SUBJECT)
    hm_temp = """# Domain Description 1.1

Interfaces 3

Interface Skull: "{0}_outer_skull_surface-low.tri"
Interface Cortex: "{0}_inner_skull_surface-low.tri"
Interface Head: "{0}_outer_skin_surface-low.tri"

Domains 4

Domain Scalp: Skull -Head
Domain Brain: -Cortex
Domain Air: Head
Domain Skull: Cortex -Skull
""".format('%s/%s/bem/watershed/%s' % (SUBJECTS_DIR, SUBJECT, SUBJECT))

    hm_geom = hm_base + '.geom'
    with open(hm_geom, 'w') as fd:
        fd.write(hm_temp)
    print ('%s written.' % (hm_geom,))

    hm_cond = hm_base + '.cond'
    with open(hm_cond, 'w') as fd:
        fd.write("""# Properties Description 1.0 (Conductivities)

Air         0.0
Scalp       1
Brain       1
Skull       0.03
""")
    print ('%s written.' % (hm_cond,))


if __name__ == '__main__':
    cmd = sys.argv[1]

    if cmd == 'gdist':
        compute_gdist_mat(*sys.argv[2:])
