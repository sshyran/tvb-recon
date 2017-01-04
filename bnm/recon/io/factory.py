# -*- coding: utf-8 -*-

import os
from bnm.recon.io.annotation import H5AnnotationIO, AnnotationIO
from bnm.recon.io.surface import GiftiSurfaceIO, FreesurferIO, H5SurfaceIO
from bnm.recon.io.volume import VolumeIO, H5VolumeIO
from bnm.recon.model.constants import GIFTI_EXTENSION, H5_EXTENSION


class IOFactory(object):
    def __get_extension(self, file_path):
        _, extension = os.path.splitext(file_path)
        return extension

    def get_surface_io(self, surface_path):
        extension = self.__get_extension(surface_path)

        if extension == GIFTI_EXTENSION:
            return GiftiSurfaceIO()
        else:
            if extension == H5_EXTENSION:
                return H5SurfaceIO()
            else:
                return FreesurferIO()

    def read_surface(self, surface_path, use_center_surface):
        surface_io = self.get_surface_io(surface_path)
        return surface_io.read(surface_path, use_center_surface)

    def write_surface(self, out_surface_path, surface):
        surface_io = self.get_surface_io(out_surface_path)
        surface_io.write(surface, out_surface_path)

    def get_volume_io(self, volume_path):
        extension = self.__get_extension(volume_path)

        if extension == H5_EXTENSION:
            return H5VolumeIO()
        else:
            return VolumeIO()

    def read_volume(self, volume_path):
        volume_io = self.get_volume_io(volume_path)
        return volume_io.read(volume_path)

    def write_volume(self, out_volume_path, volume):
        volume_io = self.get_volume_io(out_volume_path)
        volume_io.write(out_volume_path, volume)

    def get_annotation_io(self, annotation_path):
        extension = self.__get_extension(annotation_path)

        if extension == H5_EXTENSION:
            return H5AnnotationIO()
        else:
            return AnnotationIO()

    def read_annotation(self, annotation_path):
        annotation_io = self.get_annotation_io(annotation_path)
        return annotation_io.read(annotation_path)

    def write_annotation(self, out_annotation_path, annotation):
        annotation_io = self.get_annotation_io(out_annotation_path)
        annotation_io.write(out_annotation_path, annotation)
