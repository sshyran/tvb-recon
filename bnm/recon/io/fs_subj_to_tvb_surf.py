"""Generate TVB surfaces and region mapping for FreeSurfer subjects.

Usage:
    fs_subj_to_tvb_surf.py [--subject=<name>]
    fs_subj_to_tvb_surf.py (-h | --help)

Options:
    -h --help   Show this documentation
    -s          Subject name

"""

import os
import docopt
import numpy as np

from bnm.recon.io import tvb
from bnm.recon.algo.geom import merge_lh_rh
from bnm.recon.qc.parser.annotation import AnnotationIO
from bnm.recon.qc.parser.surface import FreesurferIO

if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    # prefer arg'd subject name
    if args['--subject'] is not None:
        os.environ['SUBJECT'] = args['--subject']

    subjects_dir = os.environ['SUBJECTS_DIR']
    subject = os.environ['SUBJECT']

    lh_surf_path = os.path.join(subjects_dir, subject, 'surf', 'lh.pial')
    rh_surf_path = os.path.join(subjects_dir, subject, 'surf', 'rh.pial')
    lh_annot_path = os.path.join(subjects_dir, subject, 'label', 'lh.aparc.annot')
    rh_annot_path = os.path.join(subjects_dir, subject, 'label', 'rh.aparc.annot')

    surface_io = FreesurferIO()

    # load left and right pial surfaces
    lh_surface = surface_io.read(lh_surf_path)
    rh_surface = surface_io.read(rh_surf_path)

    annot_io = AnnotationIO()

    # load left and right region mapping
    lh_annot = annot_io.read(lh_annot_path)
    rh_annot = annot_io.read(rh_annot_path)

    # merge surfaces and roi maps
    v, f, rm = merge_lh_rh(lh_surface.vertices, lh_surface.triangles, rh_surface.vertices, rh_surface.triangles,
                           lh_annot.region_mapping, rh_annot.region_mapping)

    # write out in TVB format
    np.savetxt('%s_ctx_roi_map.txt' % (subject,), rm.flat[:], '%i')
    tvb.write_surface_zip('%s_pial_surf.zip' % (subject,), v, f)
