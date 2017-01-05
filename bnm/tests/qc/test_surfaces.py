# -*- coding: utf-8 -*-
import os

import pytest
from bnm.recon.io.factory import IOFactory
from bnm.tests.base import get_data_file, get_temporary_files_path, remove_temporary_test_files
from nibabel.py3k import FileNotFoundError

TEST_FS_SUBJECT = "freesurfer_fsaverage"
TEST_MODIF_SUBJECT = "fsaverage_modified"
TEST_SURFACE_FOLDER = "surf"
io_factory = IOFactory()


def teardown_module():
    remove_temporary_test_files()


def test_parse_fs_surface():
    file_path = get_data_file(TEST_FS_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial")
    surf = io_factory.read_surface(file_path, False)
    assert len(surf.triangles) == 327680


def test_parse_fs_centered_surface():
    file_path = get_data_file(TEST_FS_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial")
    surf = io_factory.read_surface(file_path, True)
    assert len(surf.triangles) == 327680


def test_parse_centered_fs_surface():
    file_path = get_data_file(TEST_MODIF_SUBJECT, TEST_SURFACE_FOLDER, "lh-centered.pial")
    surf = io_factory.read_surface(file_path, False)
    assert surf.center_ras == [0, 0, 0]


def test_parse_not_existing_fs_surface():
    file_path = "not_existing_surface.pial"
    surface_io = io_factory.get_surface_io(file_path)
    with pytest.raises(IOError):
        surface_io.read(file_path, False)


def test_parse_not_surface():
    file_path = get_data_file(TEST_FS_SUBJECT, "label", "lh.aparc.annot")
    surface_io = io_factory.get_surface_io(file_path)
    with pytest.raises(ValueError):
        surface_io.read(file_path, False)


def test_parse_gifti_surface():
    file_path = get_data_file(TEST_MODIF_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial.gii")
    surf = io_factory.read_surface(file_path, False)
    assert len(surf.triangles) == 327680


def test_parse_gifti_centered_surface():
    file_path = get_data_file(TEST_MODIF_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial.gii")
    surf = io_factory.read_surface(file_path, True)
    assert len(surf.triangles) == 327680


def test_parse_not_existing_gifti_surface():
    file_path = "not_existing_surface.gii"
    surface_io = io_factory.get_surface_io(file_path)
    with pytest.raises(FileNotFoundError):
        surface_io.read(file_path, False)


def test_parse_h5_surface():
    h5_path = get_data_file('head2', 'SurfaceCortical.h5')
    surface = io_factory.read_surface(h5_path, False)
    assert len(surface.vertices) == 16
    assert len(surface.triangles) == 24


def test_write_fs_surface():
    file_path = get_data_file(TEST_FS_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial")
    original_surface = io_factory.read_surface(file_path, False)
    triangles_number = len(original_surface.triangles)

    output_file_path = get_temporary_files_path("lh-test.pial")
    io_factory.write_surface(output_file_path, original_surface)

    new_surface = io_factory.read_surface(output_file_path, False)
    assert triangles_number == len(new_surface.triangles) == 327680


def test_write_gifti_surface():
    file_path = get_data_file(TEST_MODIF_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial.gii")
    original_surface = io_factory.read_surface(file_path, False)
    triangles_number = len(original_surface.triangles)

    output_file_path = get_temporary_files_path("lh-test.pial.gii")
    io_factory.write_surface(output_file_path, original_surface)

    new_surface = io_factory.read_surface(output_file_path, False)
    assert triangles_number == len(new_surface.triangles) == 327680


def test_read_transformation_matrix_from_fs_metadata():
    file_path = get_data_file(TEST_FS_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial")
    surface_io = io_factory.get_surface_io(file_path)
    surf = surface_io.read(file_path, False)
    matrix = surface_io.read_transformation_matrix_from_metadata(surf.get_main_metadata())
    assert matrix.tolist() == [[-1, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0], [0, 0, 0, 1]]


def test_write_transformation_matrix_fs_metadata():
    file_path = get_data_file(TEST_FS_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial")
    surface_io = io_factory.get_surface_io(file_path)
    surf = surface_io.read(file_path, False)
    surface_io.write_transformation_matrix(surf.get_main_metadata())
    matrix = surface_io.read_transformation_matrix_from_metadata(surf.get_main_metadata())
    assert matrix.tolist() == [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]


def test_read_transformation_matrix_from_gifti_metadata():
    file_path = get_data_file(TEST_MODIF_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial.gii")
    surface_io = io_factory.get_surface_io(file_path)
    surf = surface_io.read(file_path, False)
    matrix = surface_io.read_transformation_matrix_from_metadata(surf.get_main_metadata())
    assert matrix == [[-1, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0], [0, 0, 0, 1]]


def test_write_transformation_matrix_gifti_metadata():
    file_path = get_data_file(TEST_MODIF_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial.gii")
    surface_io = io_factory.get_surface_io(file_path)
    surf = surface_io.read(file_path, False)
    surface_io.write_transformation_matrix(surf.get_main_metadata())
    matrix = surface_io.read_transformation_matrix_from_metadata(surf.get_main_metadata())
    assert matrix == [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]


def test_write_write_brain_visa_surf():
    surface_path = get_data_file(TEST_FS_SUBJECT, TEST_SURFACE_FOLDER, "lh.pial")
    out_path = get_temporary_files_path("lh.pial.tri")

    io_surface = io_factory.get_surface_io(surface_path)
    surface = io_surface.read(surface_path, False)
    io_surface.write_brain_visa_surf(out_path, surface)

    assert os.path.exists(out_path)