# -*- coding: utf-8 -*-

from ...io.factory import IOUtils
from ..base import BaseTest, get_data_file


_expected_region_names = [
    'unknown', 'bankssts', 'caudalanteriorcingulate', 'caudalmiddlefrontal',
    'corpuscallosum', 'cuneus', 'entorhinal', 'fusiform', 'inferiorparietal',
    'inferiortemporal', 'isthmuscingulate', 'lateraloccipital',
    'lateralorbitofrontal', 'lingual', 'medialorbitofrontal',
    'middletemporal', 'parahippocampal', 'paracentral', 'parsopercularis',
    'parsorbitalis', 'parstriangularis', 'pericalcarine', 'postcentral',
    'posteriorcingulate', 'precentral', 'precuneus',
    'rostralanteriorcingulate', 'rostralmiddlefrontal', 'superiorfrontal',
    'superiorparietal', 'superiortemporal', 'supramarginal', 'frontalpole',
    'temporalpole', 'transversetemporal', 'insula']


class AnnotationTests(BaseTest):

    def setUp(self):
        super().setUp()
        self.subject = 'freesurfer_fsaverage'
        self.annot_path = 'label'

    def test_parse_annotation(self):
        file_path = get_data_file(
            self.subject, self.annot_path, "lh.aparc.annot")
        annot = IOUtils.read_annotation(file_path)
        self.assertEqual(_expected_region_names, annot.region_names)

    def test_parse_not_existent_annotation(self):
        file_path = "not_existent_annotation.annot"
        annotation_io = IOUtils.annotation_io_factory(file_path)
        self.assertRaises(IOError, annotation_io.read, file_path)

    def test_parse_not_annotation(self):
        file_path = get_data_file(self.subject, "surf", "lh.pial")
        annotation_io = IOUtils.annotation_io_factory(file_path)
        self.assertRaises(ValueError, annotation_io.read, file_path)

    def test_parse_h5_annotation(self):
        h5_path = get_data_file('head2', 'RegionMapping.h5')
        annotation = IOUtils.read_annotation(h5_path)
        self.assertEqual(annotation.region_mapping.size, 16)

    def test_write_annotation(self):
        file_path = get_data_file(
            self.subject, self.annot_path, "lh.aparc.annot")
        annotation = IOUtils.read_annotation(file_path)
        out_annotation_path = self.temp_file_path("lh-test.aparc.annot")
        IOUtils.write_annotation(out_annotation_path, annotation)
        new_annotation = IOUtils.read_annotation(out_annotation_path)
        self.assertEqual(annotation.region_names, new_annotation.region_names)
