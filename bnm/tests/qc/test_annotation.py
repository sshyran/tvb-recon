# -*- coding: utf-8 -*-

import pytest
from bnm.recon.io.factory import IOFactory
from bnm.tests.base import get_data_file

TEST_SUBJECT = 'freesurfer_fsaverage'
TEST_ANNOT_FOLDER = 'label'
io_factory = IOFactory()

def test_parse_annotation():
    file_path = get_data_file(TEST_SUBJECT, TEST_ANNOT_FOLDER, "lh.aparc.annot")
    parser = io_factory.get_annotation_io(file_path)
    annot = parser.read(file_path)
    assert annot.region_names == ['unknown', 'bankssts', 'caudalanteriorcingulate', 'caudalmiddlefrontal',
                                  'corpuscallosum', 'cuneus', 'entorhinal', 'fusiform', 'inferiorparietal',
                                  'inferiortemporal', 'isthmuscingulate', 'lateraloccipital', 'lateralorbitofrontal',
                                  'lingual', 'medialorbitofrontal', 'middletemporal', 'parahippocampal', 'paracentral',
                                  'parsopercularis', 'parsorbitalis', 'parstriangularis', 'pericalcarine',
                                  'postcentral', 'posteriorcingulate', 'precentral', 'precuneus',
                                  'rostralanteriorcingulate', 'rostralmiddlefrontal', 'superiorfrontal',
                                  'superiorparietal', 'superiortemporal', 'supramarginal', 'frontalpole',
                                  'temporalpole', 'transversetemporal', 'insula']


def test_parse_not_existent_annotation():
    file_path = "not_existent_annotation.annot"
    parser = io_factory.get_annotation_io(file_path)
    with pytest.raises(IOError):
        parser.read(file_path)


def test_parse_not_annotation():
    file_path = get_data_file(TEST_SUBJECT, "surf", "lh.pial")
    parser = io_factory.get_annotation_io(file_path)
    with pytest.raises(ValueError):
        parser.read(file_path)


def test_parse_h5_annotation():
    h5_path = get_data_file('head2', 'RegionMapping.h5')
    h5_io = io_factory.get_annotation_io(h5_path)
    annotation = h5_io.read(h5_path)
    assert len(annotation.region_mapping) == 16
