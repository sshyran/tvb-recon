# -*- coding: utf-8 -*-

import pytest
from bnm.recon.qc.parser.surface import FreesurferParser, GiftiSurfaceParser
from bnm.tests.base import get_data_file, get_temporary_files_path, remove_temporary_test_files
from nibabel.filebasedimages import ImageFileError
from nibabel.py3k import FileNotFoundError

TEST_FS_SUBJECT = "freesurfer_fsaverage"
TEST_MODIF_SUBJECT = "fsaverage_modified"
TEST_SURFACE_FOLDER = "surf"


def teardown_module():
    remove_temporary_test_files()


def test_parse_fs_surface():
    parser = FreesurferParser()
    file_path = get_data_file(TEST_FS_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial")
    surf = parser.read(file_path)
    assert len(surf.triangles) == 327680


def test_parse_not_existing_fs_suraface():
    parser = FreesurferParser()
    file_path = "not_existing_surface.pial"
    with pytest.raises(IOError):
        parser.read(file_path)


def test_parse_not_fs_surface():
    parser = FreesurferParser()
    file_path = get_data_file(TEST_MODIF_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial.gii")
    with pytest.raises(ValueError):
        parser.read(file_path)


def test_parse_gifti_surface():
    parser = GiftiSurfaceParser()
    file_path = get_data_file(TEST_MODIF_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial.gii")
    surf = parser.read(file_path)
    assert len(surf.triangles) == 327680


def test_parse_not_existing_gifti_suraface():
    parser = GiftiSurfaceParser()
    file_path = "not_existing_surface.pial"
    with pytest.raises(FileNotFoundError):
        parser.read(file_path)


def test_parse_not_gifti_surface():
    parser = GiftiSurfaceParser()
    file_path = get_data_file(TEST_FS_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial")
    with pytest.raises(ImageFileError):
        parser.read(file_path)


def test_write_fs_surface():
    parser = FreesurferParser()
    file_path = get_data_file(TEST_FS_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial")
    original_surface = parser.read(file_path)
    triangles_number = len(original_surface.triangles)

    output_file_path = get_temporary_files_path("lh-test.pial")
    parser.write(original_surface, output_file_path)

    new_surface = parser.read(output_file_path)
    assert triangles_number == len(new_surface.triangles) == 327680


def test_write_gifti_surface():
    parser = GiftiSurfaceParser()
    file_path = get_data_file(TEST_MODIF_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial.gii")
    original_surface = parser.read(file_path)
    triangles_number = len(original_surface.triangles)

    output_file_path = get_temporary_files_path("lh-test.pial.gii")
    parser.write(original_surface, output_file_path)

    new_surface = parser.read(output_file_path)
    assert triangles_number == len(new_surface.triangles) == 327680


def test_read_transformation_matrix_from_fs_metadata():
    parser = FreesurferParser()
    file_path = get_data_file(TEST_FS_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial")
    surf = parser.read(file_path)
    matrix = parser.read_transformation_matrix_from_metadata(surf.get_main_metadata())
    assert matrix.tolist() == [[-1, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0], [0, 0, 0, 1]]


def test_write_transformation_matrix_fs_metadata():
    parser = FreesurferParser()
    file_path = get_data_file(TEST_FS_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial")
    surf = parser.read(file_path)
    parser.write_transformation_matrix(surf.get_main_metadata())
    matrix = parser.read_transformation_matrix_from_metadata(surf.get_main_metadata())
    assert matrix.tolist() == [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]


def test_read_transformation_matrix_from_gifti_metadata():
    parser = GiftiSurfaceParser()
    file_path = get_data_file(TEST_MODIF_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial.gii")
    surf = parser.read(file_path)
    matrix = parser.read_transformation_matrix_from_metadata(surf.get_main_metadata())
    assert matrix == [[-1, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0], [0, 0, 0, 1]]


def test_write_transformation_matrix_gifti_metadata():
    parser = GiftiSurfaceParser()
    file_path = get_data_file(TEST_MODIF_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial.gii")
    surf = parser.read(file_path)
    parser.write_transformation_matrix(surf.get_main_metadata())
    matrix = parser.read_transformation_matrix_from_metadata(surf.get_main_metadata())
    assert matrix == [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
