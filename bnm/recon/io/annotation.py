# -*- coding: utf-8 -*-

import h5py
from bnm.recon.model.annotation import Annotation
from nibabel.freesurfer.io import read_annot, write_annot


class ABCAnnotationIO(object):
    """
    This will define the behaviour needed for an annotation io.
    """

    def read(self, annotation_path):
        raise NotImplementedError()

    def write(self, out_annotation_path, annotation):
        raise NotImplementedError()


class AnnotationIO(ABCAnnotationIO):

    def read(self, annotation_path):
        region_mapping, regions_color_table, region_names = read_annot(annotation_path)
        return Annotation(region_mapping, regions_color_table, region_names)

    def write(self, out_annot_path, region_mapping, regions_color_table, region_names):
        write_annot(out_annot_path, region_mapping, regions_color_table, region_names)

    def write(self, out_annotation_path, annotation):
        write_annot(out_annotation_path, annotation.region_mapping, annotation.regions_color_table,
                    annotation.region_names)


class H5AnnotationIO(ABCAnnotationIO):

    def read(self, annotation_path):
        h5_file = h5py.File(annotation_path, 'r', libver='latest')
        region_mapping = h5_file['/data'][()]
        h5_file.close()
        return Annotation(region_mapping, [], [])
