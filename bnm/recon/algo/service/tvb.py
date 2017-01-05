# -*- coding: utf-8 -*-

import os
import numpy
from bnm.recon.io.tvb import TVBWriter


class TVBService(object):
    # Merge surfaces and roi maps. Write out in TVB format.
    def convert_fs_subj_to_tvb_surf(self, subject=None):
        subjects_dir = os.environ['SUBJECTS_DIR']

        if subject is None:
            subject = os.environ['SUBJECT']

        lh_surf_path = os.path.join(subjects_dir, subject, 'surf', 'lh.pial')
        rh_surf_path = os.path.join(subjects_dir, subject, 'surf', 'rh.pial')
        lh_annot_path = os.path.join(subjects_dir, subject, 'label', 'lh.aparc.annot')
        rh_annot_path = os.path.join(subjects_dir, subject, 'label', 'rh.aparc.annot')

        lh_surface = self.surface_io.read(lh_surf_path, False)
        rh_surface = self.surface_io.read(rh_surf_path, False)

        lh_annot = self.annotation_io.read(lh_annot_path)
        rh_annot = self.annotation_io.read(rh_annot_path)

        surface, region_mapping = self.merge_lh_rh(lh_surface, rh_surface, lh_annot.region_mapping,
                                                   rh_annot.region_mapping)

        numpy.savetxt('%s_ctx_roi_map.txt' % (subject,), region_mapping.flat[:], '%i')
        TVBWriter().write_surface_zip('%s_pial_surf.zip' % (subject,), surface)
