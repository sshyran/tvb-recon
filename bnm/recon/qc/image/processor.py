# -*- encoding: utf-8 -*-

import os
from bnm.recon.qc.image.writer import ImageWriter
from bnm.recon.qc.parser.annotation import AnnotationParser
from bnm.recon.qc.parser.generic import GenericParser
from bnm.recon.qc.parser.surface import FreesurferParser, GiftiSurfaceParser
from bnm.recon.qc.parser.volume import VolumeParser
from bnm.recon.qc.model.constants import projections


class ImageProcessor(object):
    # TODO keep a constant for .png and snapshotname?
    snapshot_name = "snapshot"
    snapshot_extension = ".png"

    def __init__(self):
        self.parser_volume = VolumeParser()
        self.generic_parser = GenericParser()
        self.annotation_parser = AnnotationParser()
        self.writer = ImageWriter()

        try:
            snapshot_count = int(os.environ['SNAPSHOT_NUMBER'])
        except ValueError:
            snapshot_count = 0
        self.snapshot_count = snapshot_count

    def generate_file_name(self, current_projection):
        file_name = self.snapshot_name + str(self.snapshot_count) + current_projection
        return file_name

    @staticmethod
    def factory_surface_parser(surface_path):
        gifti_extension = ".gii"
        filename, extension = os.path.splitext(surface_path)

        if extension == gifti_extension:
            return GiftiSurfaceParser()
        else:
            return FreesurferParser()

    def read_surface(self, surface_path):

        parser = self.factory_surface_parser(surface_path)
        return parser.read(surface_path)

    def show_single_volume(self, volume_path):

        volume = self.parser_volume.parse(volume_path)
        ras = self.generic_parser.get_ras_coordinates()

        for projection in projections:
            x_axis_coords, y_axis_coords, volume_matrix = volume.slice_volume(projection, ras)
            self.writer.write_matrix(x_axis_coords, y_axis_coords, volume_matrix, self.generate_file_name(projection))

    def overlap_2_volumes(self, background_path, overlay_path):

        background_volume = self.parser_volume.parse(background_path)
        overlay_volume = self.parser_volume.parse(overlay_path)

        ras = self.generic_parser.get_ras_coordinates()

        for projection in projections:
            x, y, background_matrix = background_volume.slice_volume(projection, ras)
            x1, y1, overlay_matrix = overlay_volume.slice_volume(projection, ras)
            self.writer.write_2_matrices(x, y, background_matrix, x1, y1, overlay_matrix,
                                         self.generate_file_name(projection))

    def overlap_3_volumes(self, background_path, overlay_1_path, overlay_2_path):

        volume_background = self.parser_volume.parse(background_path)
        volume_overlay_1 = self.parser_volume.parse(overlay_1_path)
        volume_overlay_2 = self.parser_volume.parse(overlay_2_path)

        ras = self.generic_parser.get_ras_coordinates()

        for projection in projections:
            x, y, background_matrix = volume_background.slice_volume(projection, ras)
            x1, y1, overlay_1_matrix = volume_overlay_1.slice_volume(projection, ras)
            x2, y2, overlay_2_matrix = volume_overlay_2.slice_volume(projection, ras)
            self.writer.write_3_matrices(x, y, background_matrix, x1, y1, overlay_1_matrix, x2, y2, overlay_2_matrix,
                                         self.generate_file_name(projection))

    def overlap_surface_annotation(self, surface_path, annotation):
        annotation = self.annotation_parser.parse(annotation)
        surface = self.read_surface(surface_path)
        self.writer.write_surface_with_annotation(surface, annotation, self.generate_file_name('surface_annotation'))

    def overlap_volume_surface(self, volume_background, surfaces_path):
        volume = self.parser_volume.parse(volume_background)
        # TODO review varargs processing
        surfaces = [0 for _ in xrange(len(surfaces_path))]

        for i, surf in enumerate(surfaces_path):
            surfaces[i] = self.read_surface(os.path.expandvars(surf))

        ras = self.generic_parser.get_ras_coordinates()

        for projection in projections:
            x, y, background_matrix = volume.slice_volume(projection, ras)
            clear_flag = True
            for surface in surfaces:
                x_array, y_array = surface.cut_by_plane(projection, ras)
                self.writer.write_matrix_and_surface(x, y, background_matrix, x_array, y_array, clear_flag)
                clear_flag = False
            self.writer.save_figure(self.generate_file_name(projection))

    def overlap_volume_surfaces(self, volume_background, resampled_name):
        surfaces_path = os.path.expandvars(os.environ['SURF'])
        if resampled_name != '':
            resampled_name = '.' + resampled_name
        volume = self.parser_volume.parse(volume_background)
        ras = self.generic_parser.get_ras_coordinates()
        print ras
        for i in projections:
            clear_flag = True
            x, y, background_matrix = volume.slice_volume(i, ras)
            for k in ('rh', 'lh'):
                for j in ('pial', 'white'):
                    current_surface = self.read_surface(surfaces_path + '/' + k + '.' + j + resampled_name + '.gii')
                    surf_x_array, surf_y_array = current_surface.cut_by_plane(i, ras)
                    self.writer.write_matrix_and_surfaces(x, y, background_matrix, surf_x_array, surf_y_array,
                                                          clear_flag, j)
                    clear_flag = False
            self.writer.save_figure(self.generate_file_name(i))
