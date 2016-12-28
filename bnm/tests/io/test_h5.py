# -*- coding: utf-8 -*-

from bnm.recon.io.h5 import H5IO
from bnm.tests.base import get_data_file

def test_read_surface():
    h5_path = get_data_file('head2', 'SurfaceCortical.h5')
    h5_io = H5IO()
    surface = h5_io.read_surface(h5_path)
    assert len(surface.vertices) == 16
    assert len(surface.triangles) == 24


def test_read_region_mapping():
    h5_path = get_data_file('head2', 'RegionMapping.h5')
    h5_io = H5IO()
    annotation = h5_io.read_region_mapping(h5_path)
    assert len(annotation.region_mapping) == 16

def test_read_volume():
    h5_path = get_data_file('head2', 'VolumeT1Background.h5')
    h5_io = H5IO()
    volume = h5_io.read_volume(h5_path)
    assert volume.dimensions == (6, 5, 4)