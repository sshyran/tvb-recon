# -*- coding: utf-8 -*-

import h5py
from bnm.recon.model.annotation import Annotation
from bnm.recon.model.surface import Surface
from bnm.recon.model.volume import Volume


class H5IO(object):
    def read_surface(self, h5_path):
        h5_file = h5py.File(h5_path, 'r', libver='latest')
        vertices = h5_file['/vertices'][()]
        triangles = h5_file['/triangles'][()]
        h5_file.close()
        return Surface(vertices, triangles, [], None)

    def __read_data_field(self, h5_path):
        h5_file = h5py.File(h5_path, 'r', libver='latest')
        data = h5_file['/data'][()]
        h5_file.close()
        return data

    def read_region_mapping(self, h5_path):
        region_mapping = self.__read_data_field((h5_path))
        return Annotation(region_mapping, [], [])

    def read_volume(self, h5_path):
        data = self.__read_data_field(h5_path)
        return Volume(data, [], None)