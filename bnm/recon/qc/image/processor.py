# -*- encoding: utf-8 -*-

from __future__ import print_function
import os
from .writer import ImageWriter
from ..parser import annotation, generic, surface, volume
from ..model.constants import projections


class ImageProcessor(object):
    snapshot_name = "snapshot"
    snapshot_extension = ".png"

    def __init__(self):
        self.parser_volume = volume.VolumeParser()
        self.parser_surface = surface.SurfaceParser()
        self.generic_parser = generic.GenericParser()
        self.annotation_parser = annotation.AnnotationParser()
        self.writer = ImageWriter()

        try:
            snapshot_count = int(os.environ['SNAPSHOT_NUMBER'])
        except ValueError:
            snapshot_count = 0

        self.snapshot_count = snapshot_count


    def _new_name(self, current_projection):
        file_name = self.snapshot_name + str(self.snapshot_count) + current_projection
        return file_name


    def is_surface_gifti(self, surface_path):
        gifti_extension = ".gii"
        filename, extension = os.path.splitext(surface_path)

        if(extension == gifti_extension):
            return True
        else:
            return False


    def choose_parser_for_surface(self, surface_path):

        if (self.is_surface_gifti(surface_path)):
            return self.parser_surface.parse_gifti(surface_path)
        else:
            return self.parser_surface.parse_fs(surface_path)


    def show_single_volume(self, volume_path):

        volume = self.parser_volume.parse(volume_path)
        ras = [0.68, -53.32, -19.71] #self.generic_parser.get_ras_coordinates()

        for i in projections:
            volume_matrix = volume.align(i, ras)
            self.writer.write_matrix(volume_matrix, self._new_name(i))

    def overlap_2_volumes(self, background_path, overlay_path):

        volume_background = self.parser_volume.parse(background_path)
        volume_overlay = self.parser_volume.parse(overlay_path)
        ras = [0.68, -53.32, -19.71] #self.generic_parser.get_ras_coordinates()

        for i in projections:
            background_matrix = volume_background.align(i, ras)
            overlay_matrix = volume_overlay.align(i, ras)
            self.writer.write_2_matrices(background_matrix, overlay_matrix, self._new_name(i))


    def overlap_3_volumes(self, background_path, overlay_1_path, overlay_2_path):

        volume_background = self.parser_volume.parse(background_path)
        volume_overlay_1 = self.parser_volume.parse(overlay_1_path)
        volume_overlay_2 = self.parser_volume.parse(overlay_2_path)

        ras = [0.68, -53.32, -19.71] #self.generic_parser.get_ras_coordinates()

        for i in projections:
            background_matrix = volume_background.align(i, ras)
            overlay_1_matrix = volume_overlay_1.align(i, ras)
            overlay_2_matrix = volume_overlay_2.align(i, ras)
            self.writer.write_3_matrices(background_matrix, overlay_1_matrix, overlay_2_matrix, self._new_name(i))


    def overlap_surface_annotation(self, surface_path, annotation):
        annot = self.annotation_parser.parse(annotation)
        surface = self.choose_parser_for_surface(surface_path)
        self.writer.write_surface_with_annotation(surface, annot, self._new_name('surface_annotation'))


    def overlap_volume_surface(self, volume_background, surfaces_path):
        # TODO surface contour is a little lower than it should be
        # TODO the image and the contour are reversed (compared to freeview)
        volume = self.parser_volume.parse(volume_background)
        surfaces = [0 for _ in range(len(surfaces_path))]
        for i, surf in enumerate(surfaces_path):
            surfaces[i] = self.choose_parser_for_surface(surf)
        ras = [0.68, -53.32, -19.71] #self.generic_parser.get_ras_coordinates()
        for i in projections:
            X,Y,background_matrix = volume.align(i, ras)
            clear_flag = True
            for surface in surfaces:
                x_array, y_array = surface.get_x_y_array(i, ras)
                self.writer.write_matrix_and_surface(X,Y,background_matrix, x_array, y_array, clear_flag)
                clear_flag = False
            self.writer.save_figure(self._new_name(i))


    def overlap_volume_surfaces(self, volume_background, resampled_name):
        surfaces_path = os.path.expandvars(os.environ['SURF'])
        if resampled_name != '':
            resampled_name = '.' + resampled_name
        volume = self.parser_volume.parse(volume_background)
        ras = [0.68, -53.32, -19.71] #self.generic_parser.get_ras_coordinates()
        print(ras)
        for i in projections:
            clear_flag = True
            background_matrix = volume.align(i, ras)
            for k in ('rh', 'lh'):
                for j in ('pial', 'white'):
                    current_surface = self.parser_surface.parse_gifti(surfaces_path + '/' + k + '.' + j + resampled_name + '.gii')
                    surf_x_array, surf_y_array = current_surface.get_x_y_array(i, ras)
                    self.writer.write_matrix_and_surfaces(background_matrix, surf_x_array, surf_y_array, clear_flag, j)
                    clear_flag = False
            self.writer.save_figure(self._new_name(i))

