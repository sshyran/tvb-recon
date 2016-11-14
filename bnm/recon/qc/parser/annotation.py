# -*- coding: utf-8 -*-

from nibabel.freesurfer.io import read_annot
from bnm.recon.qc.model.annotation import Annotation


class AnnotationParser(object):

    def parse(self, annotation_path):
        region_mapping, regions_color_table, region_names = read_annot(annotation_path)
        return Annotation(region_mapping, regions_color_table, region_names)
