# -*- coding: utf-8 -*-

import os
import pytest
from bnm.recon.qc.io.volume import VolumeIO
from bnm.tests.base import get_data_file, get_temporary_files_path, remove_temporary_test_files
from nibabel.filebasedimages import ImageFileError
from nibabel.py3k import FileNotFoundError

TEST_MODIF_SUBJECT = 'fsaverage_modified'
TEST_FS_SUBJECT = "freesurfer_fsaverage"
TEST_VOLUME_FOLDER = 'mri'


def teardown_module():
    remove_temporary_test_files()


def test_parse_volume():
    parser = VolumeIO()
    file_path = get_data_file(TEST_MODIF_SUBJECT, TEST_VOLUME_FOLDER, "T1.nii.gz")
    volume = parser.read(file_path)
    assert volume.dimensions == (256, 256, 256)


def test_parse_not_existing_volume():
    parser = VolumeIO()
    file_path = "not-existent-volume.nii.gz"
    with pytest.raises(FileNotFoundError):
        parser.read(file_path)


def test_parse_not_volume():
    parser = VolumeIO()
    file_path = get_data_file(TEST_FS_SUBJECT, "surf", "lh.pial")
    with pytest.raises(ImageFileError):
        parser.read(file_path)


def test_write_volume():
    parser = VolumeIO()
    in_file_path = get_data_file(TEST_MODIF_SUBJECT, TEST_VOLUME_FOLDER, "T1.nii.gz")
    volume = parser.read(in_file_path)
    out_file_path = get_temporary_files_path('T1-out.nii.gz')
    parser.write(out_file_path, volume.data, volume.affine_matrix, volume.header)
    assert os.path.exists(out_file_path)
