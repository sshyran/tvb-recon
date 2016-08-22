import os
import sys
import glob
import numpy as np
import nibabel
import scipy.io
import scipy.cluster
import warnings
import nibabel.freesurfer


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
    v, f = nibabel.freesurfer.read_geometry(fs_surf)
    write_brain_visa_surf(fs_surf + '.tri', v, f)


def annot_path(hemi, annot_name):
    annot_fname = '%s.%s.annot' % (hemi, annot_name)
    annot_path = os.path.join(SUBJECTS_DIR, SUBJECT, 'label', annot_fname)
    return annot_path


def read_annot(hemi, annot_name):
    return nibabel.freesurfer.read_annot(annot_path(hemi, annot_name))


def write_annot(hemi, annot_name, labels, ctab, names):
    return nibabel.freesurfer.write_annot(annot_path(hemi, annot_name), labels, ctab, names)


def read_surf(hemi, name):
    surf_fname = '%s.%s' % (hemi, name)
    surf_path = os.path.join(SUBJECTS_DIR, SUBJECT, 'surf', surf_fname)
    return nibabel.freesurfer.read_geometry(surf_path)


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


def make_subparc(v, f, annot, roi_names, trg_area=100.0):

    # build vertex -> face list map
    vfm = [set([]) for _ in v]
    for i, face in enumerate(f):
        for j in face:
            vfm[j].add(i)
    vfm = np.array(vfm)

    # make new annotation
    new_annot = annot.copy()
    new_names = [] # such that new_names[new_annot[i]] is correct name for i'th vertex
    next_aval = 1
    for i in np.unique(annot):
        name = roi_names[i]
        mask = annot == i

        # "unknown", just skip
        if i == -1:
            new_annot[mask] = 0
            new_names.append(name)
            continue

        # indices of faces in ROI
        rfi = set([])
        for face_set in vfm[mask]:
            rfi.update(face_set)
        rfi = np.array(list(rfi))

        # empty roi
        if rfi.size == 0:
            continue

        # compute area of faces in roi
        tri_xyz = v[f[rfi]]
        i, j, k = np.transpose(tri_xyz, (1, 0, 2))
        ij = j - i
        ik = k - i
        roi_area = np.sum(np.sqrt(np.sum(np.cross(ij, ik)**2, axis=1)) / 2.0)#}}}

        # choose k for desired roi area
        k = int(roi_area / trg_area)

        # cluster centered vertices
        v_roi = v[mask]
        _, i_lab = scipy.cluster.vq.kmeans2(v_roi - v_roi.mean(axis=0), k)

        # update annot
        new_annot[mask] = next_aval + i_lab
        next_aval += k
        new_names += ['%s-%d' % (name.decode('ascii'), j) for j in range(k)]

    # create random colored ctab
    new_ctab = np.random.randint(255, size=(len(new_names), 5))
    r, g, b, _, _ = new_ctab.T
    new_ctab[:, 3] = 0
    new_ctab[:, 4] = r + 256 * g + 256 * 256 * b # fs magic values

    return new_annot, new_ctab, new_names


def subparc_files(hemi, parc_name, out_parc_name, trg_area):
    trg_area = float(trg_area)
    v, f = read_surf(hemi, 'sphere')
    lab, ctab, names = read_annot(hemi, parc_name)
    new_lab, new_ctab, new_names = make_subparc(v, f, lab, names)
    write_annot(hemi, out_parc_name, new_lab, new_ctab, new_names)


if __name__ == '__main__':
    cmd = sys.argv[1]

    if cmd == 'gdist':
        compute_gdist_mat(*sys.argv[2:])
    if cmd == 'subparc':
        subparc_files(*sys.argv[2:])

