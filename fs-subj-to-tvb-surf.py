
import os
import os.path
import numpy as np
import nibabel.freesurfer
import zipfile


try:
    from cStringIO import StringIO
except ImportError: # Py 3
    from io import BytesIO as StringIO


def get_default_base_path():
    base_path = os.path.join(
            os.environ['SUBJECTS_DIR'],
            os.environ['SUBJECT']
        )
    return base_path


def load_surf(base_path, name):
    return nibabel.freesurfer.read_geometry(
            os.path.join(base_path, 'surf/' + name))


def load_roi_map(base_path, name):
    return nibabel.freesurfer.read_annot(
            os.path.join(base_path, 'label/' + name))[0]


def np_save_strio(arr, fmt):
    sio = StringIO()
    np.savetxt(sio, arr, fmt)
    return sio


def write_surface_zip(zip_fname, v, f):
    sv = np_save_strio(v, '%f')
    sf = np_save_strio(f, '%d')
    szf = StringIO()
    zf = zipfile.ZipFile(szf, 'w')
    zf.writestr('vertices.txt', sv.getvalue())
    zf.writestr('triangles.txt', sf.getvalue())
    zf.close()
    with open(zip_fname, 'wb') as fd:
        fd.write(szf.getvalue())


if __name__ == '__main__':
    path = get_default_base_path()

    # load left and right pial surfaces
    lv, lf = load_surf(path, 'lh.pial')
    rv, rf = load_surf(path, 'rh.pial')

    # load left and right region mapping
    lrm = load_roi_map(path, 'lh.aparc.annot')
    rrm = load_roi_map(path, 'rh.aparc.annot')

    # merge surfaces and roi maps
    v = np.r_[lv, rv]
    f = np.r_[lf, rf + lf.max()]
    rm = np.r_[lrm, rrm + lf.max()]

    # write out in TVB format
    np.savetxt('region_mapping.txt', rm.flat[:], '%i')
    write_surface_zip('surface.zip', v, f)
