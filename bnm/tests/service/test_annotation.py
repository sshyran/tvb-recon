# -*- coding: utf-8 -*-

import os
from collections import OrderedDict
from bnm.tests.base import get_data_file, get_temporary_files_path, remove_temporary_test_files

os.environ["FREESURFER_HOME"] = ""
from bnm.recon.algo.service.annotation import AnnotationService

color_lut_name = "colorLUT.txt"


def teardown_module():
    remove_temporary_test_files()


def test_read_lut_by_name():
    lut_path = get_data_file(color_lut_name)
    service = AnnotationService()
    labels, names, colors = service.read_lut(lut_path=lut_path, key_mode='name')
    assert isinstance(labels, OrderedDict)
    assert isinstance(colors, OrderedDict)
    assert names == labels.keys() == colors.keys() == ['Unknown', 'Left-Cerebral-Exterior',
                                                       'Left-Cerebral-White-Matter', 'Left-Thalamus-Proper',
                                                       'Left-Caudate', 'ctx-lh-unknown', 'ctx-lh-bankssts',
                                                       'ctx-rh-unknown', 'ctx-rh-bankssts']


def test_read_lut_by_label():
    lut_path = get_data_file(color_lut_name)
    service = AnnotationService()
    labels, names, colors = service.read_lut(lut_path=lut_path, key_mode='label')
    assert isinstance(names, OrderedDict)
    assert isinstance(colors, OrderedDict)
    assert labels == names.keys() == colors.keys() == [0, 1, 2, 10, 11, 1000, 1001, 2000, 2001]


def test_rgb_to_fs_magic_number():
    service = AnnotationService()
    fs_magic_number = service.rgb_to_fs_magic_number([100, 100, 100])
    assert fs_magic_number == 6579300


def test_annot_to_lut():
    service = AnnotationService()
    lut_path = get_temporary_files_path('colorLUT-temp.txt')
    service.annot_to_lut(get_data_file('freesurfer_fsaverage', 'label', 'lh.aparc.annot'), lut_path=lut_path)
    assert os.path.exists(lut_path)


def test_lut_to_annot_names_ctab():
    service = AnnotationService()
    lut_path = get_data_file(color_lut_name)
    names1, ctab1 = service.lut_to_annot_names_ctab(lut_path=lut_path)
    names2, ctab2 = service.lut_to_annot_names_ctab(lut_path=lut_path, labels=[0, 1, 2])
    names3, ctab3 = service.lut_to_annot_names_ctab(lut_path=lut_path, labels="0 1 2")
    assert names1[0] == names2[0] == names3[0] == 'Unknown'
    assert all(ctab1[0]) == all(ctab2[0]) == all(ctab3[0]) == all(
        [0, 0, 0, 0, service.rgb_to_fs_magic_number([0, 0, 0])])


def test_annot_names_to_labels():
    service = AnnotationService()
    lut_path = get_data_file(color_lut_name)
    labels = service.annot_names_to_labels(['Unknown'], lut_path=lut_path)
    assert labels[0] == 0
    labels = service.annot_names_to_labels(['unknown'], add_string='ctx-lh-', lut_path=lut_path)
    assert labels[0] == 1000
    labels = service.annot_names_to_labels(['bankssts'], add_string='ctx-rh-', lut_path=lut_path)
    assert labels[0] == 2001


def test_read_input_labels():
    service = AnnotationService()
    labels1, length_labels1 = service.read_input_labels(labels='0 1', ctx='rh')
    labels2, length_labels2 = service.read_input_labels(labels='0 1', ctx='lh')
    assert length_labels1 == length_labels2 == 38
    labels, length_labels = service.read_input_labels()
    assert labels == []
    assert length_labels == 0
