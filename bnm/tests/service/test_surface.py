# -*- coding: utf-8 -*-

import os
import numpy
from bnm.recon.algo.service.surface import SurfaceService
from bnm.recon.io.annotation import AnnotationIO
from bnm.recon.io.factory import IOUtils
from bnm.recon.io.surface import FreesurferIO, H5SurfaceIO
from bnm.recon.model.surface import Surface
from bnm.tests.base import get_data_file, get_temporary_files_path, remove_temporary_test_files

service = SurfaceService()


def teardown_module():
    remove_temporary_test_files()


def test_tri_area():
    h5_path = get_data_file('head2', 'SurfaceCortical.h5')
    h5_io = H5SurfaceIO()
    surface = h5_io.read(h5_path)
    rfi = numpy.array([0])
    area = service.tri_area(surface.vertices[surface.triangles[rfi]])
    assert area[0] == 5000


def test_merge_lh_rh():
    h5_surface_path = get_data_file("head2", "SurfaceCortical.h5")
    h5_surface = IOUtils.read_surface(h5_surface_path, False)

    lh_surface = Surface(h5_surface.vertices[:len(h5_surface.vertices) / 2],
                         h5_surface.triangles[:len(h5_surface.triangles) / 2], [], None)
    rh_surface = Surface(h5_surface.vertices[len(h5_surface.vertices) / 2:],
                         h5_surface.triangles[len(h5_surface.triangles) / 2:], [], None)

    h5_region_mapping_path = get_data_file("head2", "RegionMapping.h5")
    annotation = IOUtils.read_annotation(h5_region_mapping_path)

    lh_region_mapping = annotation.region_mapping[:len(annotation.region_mapping) / 2]
    rh_region_mapping = annotation.region_mapping[len(annotation.region_mapping) / 2:]

    out_surface, out_region_mapping = service.merge_lh_rh(lh_surface, rh_surface, lh_region_mapping, rh_region_mapping)

    assert len(out_surface.vertices) == len(lh_surface.vertices) + len(rh_surface.vertices)
    assert len(out_surface.triangles) == len(lh_surface.triangles) + len(rh_surface.triangles)
    assert len(out_region_mapping) == len(lh_region_mapping) + len(rh_region_mapping)


def test_extract_subsurf():
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

    annotation = IOUtils.read_annotation(out_annot_path)
    assert numpy.array_equal(annotation.regions_color_table,
                             [[0, 118, 14, 255, 947712], [122, 186, 220, 255, 14465658]])


def test_vertex_connectivity():
    surf_path = get_data_file("head2", 'SurfaceCortical.h5')
    surface_parser = H5SurfaceIO()
    surface = surface_parser.read(surf_path)
    conn = service.vertex_connectivity(surface.vertices, surface.triangles)
    assert conn.shape == (16, 16)
    assert conn[0, 1] == 1
    assert conn[0, 10] == 0


def test_vertex_connectivity_2():
    surf_path = get_data_file("head2", 'SurfaceCortical.h5')
    surface_parser = H5SurfaceIO()
    surface = surface_parser.read(surf_path)
    conn = service.vertex_connectivity(surface.vertices, surface.triangles, metric='euclidean')
    assert conn.shape == (16, 16)
    assert conn[0, 1] == 100
    assert conn[0, 10] == 0
