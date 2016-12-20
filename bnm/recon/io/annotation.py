# -*- coding: utf-8 -*-

from nibabel.freesurfer.io import read_annot, write_annot
from bnm.recon.qc.model.annotation import Annotation


class AnnotationIO(object):

    def read(self, annotation_path):
        region_mapping, regions_color_table, region_names = read_annot(annotation_path)
        return Annotation(region_mapping, regions_color_table, region_names)

    def write(self, out_annot_path, region_mapping, regions_color_table, region_names):
        write_annot(out_annot_path, region_mapping, regions_color_table, region_names)

    def write(self, out_annotation_path, annotation):
        write_annot(out_annotation_path, annotation.region_mapping, annotation.regions_color_table,
                    annotation.region_names)
