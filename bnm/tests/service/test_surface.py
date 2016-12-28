# -*- coding: utf-8 -*-

import os
import numpy
from bnm.recon.algo.service.surface import SurfaceService
from bnm.recon.io.annotation import AnnotationIO
from bnm.recon.io.h5 import H5IO
from bnm.recon.io.surface import FreesurferIO
from bnm.tests.base import get_data_file, get_temporary_files_path, remove_temporary_test_files


def teardown_module():
    remove_temporary_test_files()

def test_tri_area():
    service = SurfaceService()
    h5_path = get_data_file('head2', 'SurfaceCortical.h5')
    h5_io = H5IO()
    surface = h5_io.read_surface(h5_path)
    rfi = numpy.array([0])
    area = service.tri_area(surface.vertices[surface.triangles[rfi]])
    assert area[0] == 5000

def test_extract_subsurf():
    service = SurfaceService()
    surface_parser = FreesurferIO()
    annot_parser = AnnotationIO()
    surface_file = get_data_file("freesurfer_fsaverage", "surf", "lh.pial")
    annot_file = get_data_file("freesurfer_fsaverage", "label", "lh.aparc.annot")
    surface = surface_parser.read(surface_file, False)
    verts = surface.vertices
    faces = surface.triangles
    annot = annot_parser.read(annot_file)
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

    surface_parser = FreesurferIO()
    surface = surface_parser.read(out_surf_path, False)
    assert len(surface.vertices) == 5714
    assert len(surface.triangles) == 11420


def test_vertex_connectivity():
    service = SurfaceService()
    surf_path = get_data_file("head2", 'SurfaceCortical.h5')
    surface_parser = H5IO()
    surface = surface_parser.read_surface(surf_path)
    conn = service.vertex_connectivity(surface.vertices, surface.triangles)
    assert conn.shape == (16, 16)
    assert conn[0, 1] == 1
    assert conn[0, 10] == 0

def test_vertex_connectivity_2():
    service = SurfaceService()
    surf_path = get_data_file("head2", 'SurfaceCortical.h5')
    surface_parser = H5IO()
    surface = surface_parser.read_surface(surf_path)
    conn = service.vertex_connectivity(surface.vertices, surface.triangles, metric='euclidean')
    assert conn.shape == (16, 16)
    assert conn[0, 1] == 100
    assert conn[0, 10] == 0