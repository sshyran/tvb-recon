# -*- coding: utf-8 -*-

import os
import pytest
from bnm.recon.io.factory import IOFactory
from bnm.tests.base import get_data_file, get_temporary_files_path, remove_temporary_test_files
from nibabel.filebasedimages import ImageFileError
from nibabel.py3k import FileNotFoundError

TEST_MODIF_SUBJECT = 'fsaverage_modified'
TEST_FS_SUBJECT = "freesurfer_fsaverage"
TEST_VOLUME_FOLDER = 'mri'
io_factory = IOFactory()


def teardown_module():
    remove_temporary_test_files()


def test_parse_volume():
    file_path = get_data_file(TEST_MODIF_SUBJECT, TEST_VOLUME_FOLDER, "T1.nii.gz")
    volume = io_factory.read_volume(file_path)
    assert volume.dimensions == (256, 256, 256)


def test_parse_not_existing_volume():
    file_path = "not-existent-volume.nii.gz"
    volume_io = io_factory.get_volume_io(file_path)
    with pytest.raises(FileNotFoundError):
        volume_io.read(file_path)


def test_parse_not_volume():
    file_path = get_data_file(TEST_FS_SUBJECT, "surf", "lh.pial")
    volume_io = io_factory.get_volume_io(file_path)
    with pytest.raises(ImageFileError):
        volume_io.read(file_path)


def test_parse_h5_volume():
    h5_path = get_data_file('head2', 'VolumeT1Background.h5')
    volume = io_factory.read_volume(h5_path)
    assert volume.dimensions == (6, 5, 4)


def test_write_volume():
    in_file_path = get_data_file(TEST_MODIF_SUBJECT, TEST_VOLUME_FOLDER, "T1.nii.gz")
    volume = io_factory.read_volume(in_file_path)
    out_file_path = get_temporary_files_path('T1-out.nii.gz')
    io_factory.write_volume(out_file_path, volume)
    assert os.path.exists(out_file_path)
