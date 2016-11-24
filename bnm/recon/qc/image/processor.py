# -*- encoding: utf-8 -*-

import os
from bnm.recon.qc.image.writer import ImageWriter
from bnm.recon.qc.parser.annotation import AnnotationParser
from bnm.recon.qc.parser.generic import GenericParser
from bnm.recon.qc.parser.surface import FreesurferParser, GiftiSurfaceParser
from bnm.recon.qc.parser.volume import VolumeParser
from bnm.recon.qc.model.constants import PROJECTIONS, SNAPSHOT_NAME, GIFTI_EXTENSION


class ImageProcessor(object):
    def __init__(self, snapshots_directory, snapshot_count=0):
        self.parser_volume = VolumeParser()
        self.generic_parser = GenericParser()
        self.annotation_parser = AnnotationParser()
        self.writer = ImageWriter(snapshots_directory)
        self.snapshot_count = snapshot_count

    def generate_file_name(self, current_projection):
        file_name = SNAPSHOT_NAME + str(self.snapshot_count) + current_projection
        return file_name

    @staticmethod
    def factory_surface_parser(surface_path):
        filename, extension = os.path.splitext(surface_path)

        if extension == GIFTI_EXTENSION:
            return GiftiSurfaceParser()
        else:
            return FreesurferParser()

    def read_surface(self, surface_path, use_center_surface):

        parser = self.factory_surface_parser(surface_path)
        return parser.read(surface_path, use_center_surface)

    def show_single_volume(self, volume_path):

        volume = self.parser_volume.parse(volume_path)
        ras = self.generic_parser.get_ras_coordinates()

        for projection in PROJECTIONS:
            x_axis_coords, y_axis_coords, volume_matrix = volume.slice_volume(projection, ras)
            self.writer.write_matrix(x_axis_coords, y_axis_coords, volume_matrix, self.generate_file_name(projection))

    def overlap_2_volumes(self, background_path, overlay_path):

        background_volume = self.parser_volume.parse(background_path)
        overlay_volume = self.parser_volume.parse(overlay_path)

        ras = self.generic_parser.get_ras_coordinates()

        for projection in PROJECTIONS:
            x, y, background_matrix = background_volume.slice_volume(projection, ras)
            x1, y1, overlay_matrix = overlay_volume.slice_volume(projection, ras)
            self.writer.write_2_matrices(x, y, background_matrix, x1, y1, overlay_matrix,
                                         self.generate_file_name(projection))

    def overlap_3_volumes(self, background_path, overlay_1_path, overlay_2_path):

        volume_background = self.parser_volume.parse(background_path)
        volume_overlay_1 = self.parser_volume.parse(overlay_1_path)
        volume_overlay_2 = self.parser_volume.parse(overlay_2_path)

        ras = self.generic_parser.get_ras_coordinates()

        for projection in PROJECTIONS:
            x, y, background_matrix = volume_background.slice_volume(projection, ras)
            x1, y1, overlay_1_matrix = volume_overlay_1.slice_volume(projection, ras)
            x2, y2, overlay_2_matrix = volume_overlay_2.slice_volume(projection, ras)
            self.writer.write_3_matrices(x, y, background_matrix, x1, y1, overlay_1_matrix, x2, y2, overlay_2_matrix,
                                         self.generate_file_name(projection))

    def overlap_surface_annotation(self, surface_path, annotation):
        annotation = self.annotation_parser.parse(annotation)
        surface = self.read_surface(surface_path, False)
        self.writer.write_surface_with_annotation(surface, annotation, self.generate_file_name('surface_annotation'))

    def overlap_volume_surfaces(self, volume_background, surfaces_path, use_center_surface):
        volume = self.parser_volume.parse(volume_background)
        ras = self.generic_parser.get_ras_coordinates()
        surfaces = [self.read_surface(os.path.expandvars(surface), use_center_surface) for surface in surfaces_path]

        for projection in PROJECTIONS:
            x, y, background_matrix = volume.slice_volume(projection, ras)
            clear_flag = True
            for surface_index, surface in enumerate(surfaces):
                surf_x_array, surf_y_array = surface.cut_by_plane(projection, ras)
                self.writer.write_matrix_and_surfaces(x, y, background_matrix, surf_x_array, surf_y_array,
                                                      surface_index, clear_flag)
                clear_flag = False
            self.writer.save_figure(self.generate_file_name(projection))
