# -*- coding: utf-8 -*-
import os

import numpy
import pytest
from bnm.tests.base import get_data_file, get_temporary_files_path, remove_temporary_test_files
from bnm.recon.qc.parser.surface import FreesurferParser
from bnm.recon.qc.parser.annotation import AnnotationParser
from bnm.recon.algo.service.surface import SurfaceService


def teardown_module():
    remove_temporary_test_files()


def test_extract_subsurf():
    service = SurfaceService()
    surface_parser = FreesurferParser()
    annot_parser = AnnotationParser()
    surface_file = get_data_file("freesurfer_fsaverage", "surf", "lh.pial")
    annot_file = get_data_file("freesurfer_fsaverage", "label", "lh.aparc.annot")
    surface = surface_parser.read(surface_file, False)
    verts = surface.vertices
    faces = surface.triangles
    annot = annot_parser.parse(annot_file)
    labels = annot.region_mapping
    verts_mask = labels == 7
    subsurf_verts, subsurf_faces = service.extract_subsurf(verts, faces, verts_mask)
    assert subsurf_verts.all() == verts[numpy.where(verts_mask)].all()


def test_aseg_surf_conc_annot():
    service = SurfaceService()
    surf_path = "data/aseg"
    out_surf_path = get_temporary_files_path("out_aseg")
    out_annot_path = get_temporary_files_path("out_annot")
    labels = "10 11"
    colorLUT = get_data_file("colorLUT.txt")
    service.aseg_surf_conc_annot(surf_path, out_surf_path, out_annot_path, labels, colorLUT)
    assert os.path.exists(out_surf_path)
    assert os.path.exists(out_annot_path)

    surface_parser = FreesurferParser()
    surface = surface_parser.read(out_surf_path, False)
    assert len(surface.vertices) == 5714
    assert len(surface.triangles) == 11420


@pytest.mark.skip("Because we need simpler test data")
def test_vertex_connectivity():
    service = SurfaceService()
    surf_path = get_data_file("freesurfer_fsaverage", "surf", "lh.pial")
    surface_parser = FreesurferParser()
    surface = surface_parser.read(surf_path, False)
    service.vertex_connectivity(surface.vertices, surface.triangles)