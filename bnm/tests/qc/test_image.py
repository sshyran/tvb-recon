# -*- coding: utf-8 -*-

import os
import numpy
from bnm.recon.io.factory import IOUtils
from bnm.recon.model.constants import SNAPSHOTS_DIRECTORY, SNAPSHOT_NAME, AXIAL
from bnm.recon.qc.image.processor import ImageProcessor
from bnm.recon.qc.image.writer import ImageWriter
from bnm.tests.base import get_data_file

SNAPSHOT_NUMBER = 0
TEST_SUBJECT_MODIF = "fsaverage_modified"
TEST_SUBJECT = "freesurfer_fsaverage"
TEST_HEAD2 = "head2"
TEST_MRI_FOLDER = "mri"
TEST_SURF_FOLDER = "surf"
TEST_ANNOT_FOLDER = "label"
TEST_T1 = "T1.nii.gz"
TEST_BRAIN = "brain.nii.gz"
TEST_SURF = "lh.pial"
TEST_GIFTI_SURF = "lh.pial.gii"
TEST_ANNOT = "lh.aparc.annot"
processor = ImageProcessor(SNAPSHOTS_DIRECTORY, SNAPSHOT_NUMBER)


def setup_module():
    os.environ["MRI"] = os.path.join("data", TEST_SUBJECT_MODIF, TEST_MRI_FOLDER)
    os.environ["T1_RAS"] = TEST_T1
    os.environ["SUBJ_DIR"] = os.path.join("data", TEST_SUBJECT_MODIF)


def teardown_module():
    for file_path in os.listdir(SNAPSHOTS_DIRECTORY):
        os.remove(os.path.join(SNAPSHOTS_DIRECTORY, file_path))
    os.rmdir(SNAPSHOTS_DIRECTORY)


def test_show_single_volume():
    volume_path = get_data_file(TEST_SUBJECT_MODIF, TEST_MRI_FOLDER, TEST_T1)
    processor.show_single_volume(volume_path, False)
    resulted_file_name = processor.generate_file_name(AXIAL, SNAPSHOT_NAME)
    assert os.path.exists(processor.writer.get_path(resulted_file_name))


def test_overlap_2_volumes():
    background_path = get_data_file(TEST_SUBJECT_MODIF, TEST_MRI_FOLDER, TEST_T1)
    overlay_path = get_data_file(TEST_SUBJECT_MODIF, TEST_MRI_FOLDER, TEST_BRAIN)
    processor.overlap_2_volumes(background_path, overlay_path, False)
    resulted_file_name = processor.generate_file_name(AXIAL, SNAPSHOT_NAME)
    assert os.path.exists(processor.writer.get_path(resulted_file_name))


def test_overlap_3_volumes():
    background_path = get_data_file(TEST_SUBJECT_MODIF, TEST_MRI_FOLDER, TEST_T1)
    overlay1_path = get_data_file(TEST_SUBJECT_MODIF, TEST_MRI_FOLDER, TEST_BRAIN)
    overlay2_path = get_data_file(TEST_SUBJECT_MODIF, TEST_MRI_FOLDER, TEST_BRAIN)
    processor.overlap_3_volumes(background_path, overlay1_path, overlay2_path, False)
    resulted_file_name = processor.generate_file_name(AXIAL, SNAPSHOT_NAME)
    assert os.path.exists(processor.writer.get_path(resulted_file_name))


def test_overlap_surface_annotation():
    writer = ImageWriter(SNAPSHOTS_DIRECTORY)
    surface_path = get_data_file(TEST_HEAD2, "SurfaceCortical.h5")
    surface = IOUtils.read_surface(surface_path, False)
    annot_path = get_data_file(TEST_HEAD2, "RegionMapping.h5")
    annot = IOUtils.read_annotation(annot_path)
    annot.region_names = ['reg1', 'reg2']
    annot.regions_color_table = numpy.array([[200, 200, 200, 255, 30567], [100, 150, 200, 255, 30568]])
    resulted_file_name = processor.generate_file_name('surface_annotation', SNAPSHOT_NAME)
    writer.write_surface_with_annotation(surface, annot, resulted_file_name)
    assert os.path.exists(writer.get_path(resulted_file_name + str(0)))


def test_overlap_volume_surface():
    volume_path = get_data_file(TEST_SUBJECT_MODIF, TEST_MRI_FOLDER, TEST_T1)
    surface_path = get_data_file(TEST_SUBJECT, TEST_SURF_FOLDER, TEST_SURF)
    processor.overlap_volume_surfaces(volume_path, [surface_path], False, False)
    resulted_file_name = processor.generate_file_name(AXIAL, SNAPSHOT_NAME)
    assert os.path.exists(processor.writer.get_path(resulted_file_name))


def test_overlap_volume_centered_surface():
    volume_path = get_data_file(TEST_SUBJECT_MODIF, TEST_MRI_FOLDER, TEST_T1)
    surface_path = get_data_file(TEST_SUBJECT, TEST_SURF_FOLDER, TEST_SURF)
    processor.overlap_volume_surfaces(volume_path, [surface_path], True, False)
    resulted_file_name = processor.generate_file_name(AXIAL, SNAPSHOT_NAME)
    assert os.path.exists(processor.writer.get_path(resulted_file_name))


def test_overlap_volume_gifti_surface():
    volume_path = get_data_file(TEST_SUBJECT_MODIF, TEST_MRI_FOLDER, TEST_T1)
    surface_path = get_data_file(TEST_SUBJECT_MODIF, TEST_SURF_FOLDER, TEST_GIFTI_SURF)
    processor.overlap_volume_surfaces(volume_path, [surface_path], False, False)
    resulted_file_name = processor.generate_file_name(AXIAL, SNAPSHOT_NAME)
    assert os.path.exists(processor.writer.get_path(resulted_file_name))


def test_show_aparc_aseg_with_new_values():
    volume_path = get_data_file(TEST_SUBJECT_MODIF, TEST_MRI_FOLDER, "aparc+aseg.nii.gz")
    conn_measure_path = get_data_file("connectivity_measure.txt")
    processor.show_aparc_aseg_with_new_values(volume_path, conn_measure_path, '', False)
    resulted_file_name = processor.generate_file_name(AXIAL, SNAPSHOT_NAME)
    assert os.path.exists(processor.writer.get_path(resulted_file_name))
