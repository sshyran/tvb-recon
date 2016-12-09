
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

from bnm.recon.io import fs, tvb
from bnm.recon.algo.geom import merge_lh_rh

if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    # prefer arg'd subject name
    if args['--subject'] is not None:
        os.environ['SUBJECT'] = args['--subject']

    fsio = fs.FreeSurferIO.from_env_vars()

    # load left and right pial surfaces
    lv, lf = fsio.read_surf('lh', 'pial')
    rv, rf = fsio.read_surf('rh', 'pial')

    # load left and right region mapping
    lrm = fsio.read_annot('lh', 'aparc')[0]
    rrm = fsio.read_annot('rh', 'aparc')[0]

    # merge surfaces and roi maps
    v, f, rm = merge_lh_rh(lv, lf, rv, rf, lrm, rrm)

    # write out in TVB format
    np.savetxt('%s_ctx_roi_map.txt' % (fsio.subject, ), rm.flat[:], '%i')
    tvb.write_surface_zip('%s_pial_surf.zip' % (fsio.subject, ), v, f)
