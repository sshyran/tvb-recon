# -*- encoding: utf-8 -*-

import os
import numpy
from bnm.recon.io.annotation import AnnotationIO
from bnm.recon.io.generic import GenericIO
from bnm.recon.io.surface import FreesurferIO, GiftiSurfaceIO
from bnm.recon.io.volume import VolumeIO
from bnm.recon.logger import get_logger
from bnm.recon.model.constants import PROJECTIONS, SNAPSHOT_NAME, GIFTI_EXTENSION, T1_RAS_VOLUME, MRI_DIRECTORY, \
    FS_TO_CONN_INDICES_MAPPING_PATH
from bnm.recon.qc.image.writer import ImageWriter


class ImageProcessor(object):
    def __init__(self, snapshots_directory, snapshot_count=0):
        self.volume_io = VolumeIO()
        self.generic_io = GenericIO()
        self.annotation_io = AnnotationIO()
        self.writer = ImageWriter(snapshots_directory)
        self.snapshot_count = snapshot_count
        self.logger = get_logger(__name__)

    def generate_file_name(self, current_projection, snapshot_name=SNAPSHOT_NAME):
        file_name = snapshot_name + "_" + current_projection + "_" + str(self.snapshot_count)
        return file_name

    def read_t1_affine_matrix(self):
        t1_volume = self.volume_io.read(os.path.join(os.environ[MRI_DIRECTORY], os.environ[T1_RAS_VOLUME]))
        return t1_volume.affine_matrix

    @staticmethod
    def factory_surface_io(surface_path):
        filename, extension = os.path.splitext(surface_path)

        if extension == GIFTI_EXTENSION:
            return GiftiSurfaceIO()
        else:
            return FreesurferIO()

    def read_surface(self, surface_path, use_center_surface):

        surface_io = self.factory_surface_io(surface_path)
        return surface_io.read(surface_path, use_center_surface)

    def show_single_volume(self, volume_path, use_cc_point, snapshot_name=SNAPSHOT_NAME):

        volume = self.volume_io.read(volume_path)

        if use_cc_point:
            ras = self.generic_io.get_ras_coordinates(self.read_t1_affine_matrix())
        else:
            ras = volume.get_center_point()

        for projection in PROJECTIONS:
            try:
                x_axis_coords, y_axis_coords, volume_matrix = volume.slice_volume(projection, ras)
            except IndexError:
                new_ras = volume.get_center_point()
                x_axis_coords, y_axis_coords, volume_matrix = volume.slice_volume(projection, new_ras)
                self.logger.info("The volume center point has been used for %s snapshot of %s.", projection,
                                 volume_path)

            self.writer.write_matrix(x_axis_coords, y_axis_coords, volume_matrix,
                                     self.generate_file_name(projection, snapshot_name))

    def overlap_2_volumes(self, background_path, overlay_path, use_cc_point, snapshot_name=SNAPSHOT_NAME):

        background_volume = self.volume_io.read(background_path)
        overlay_volume = self.volume_io.read(overlay_path)

        if use_cc_point:
            ras = self.generic_io.get_ras_coordinates(self.read_t1_affine_matrix())
        else:
            ras = background_volume.get_center_point()

        for projection in PROJECTIONS:
            try:
                x, y, background_matrix = background_volume.slice_volume(projection, ras)
                x1, y1, overlay_matrix = overlay_volume.slice_volume(projection, ras)
            except IndexError:
                new_ras = background_volume.get_center_point()
                x, y, background_matrix = background_volume.slice_volume(projection, new_ras)
                x1, y1, overlay_matrix = overlay_volume.slice_volume(projection, new_ras)
                self.logger.info("The volume center point has been used for %s snapshot of %s and %s.", projection,
                                 background_path, overlay_path)

            self.writer.write_2_matrices(x, y, background_matrix, x1, y1, overlay_matrix,
                                         self.generate_file_name(projection, snapshot_name))

    def overlap_3_volumes(self, background_path, overlay_1_path, overlay_2_path, use_cc_point,
                          snapshot_name=SNAPSHOT_NAME):

        volume_background = self.volume_io.read(background_path)
        volume_overlay_1 = self.volume_io.read(overlay_1_path)
        volume_overlay_2 = self.volume_io.read(overlay_2_path)

        if use_cc_point:
            ras = self.generic_io.get_ras_coordinates(self.read_t1_affine_matrix())
        else:
            ras = volume_background.get_center_point()

        for projection in PROJECTIONS:
            try:
                x, y, background_matrix = volume_background.slice_volume(projection, ras)
                x1, y1, overlay_1_matrix = volume_overlay_1.slice_volume(projection, ras)
                x2, y2, overlay_2_matrix = volume_overlay_2.slice_volume(projection, ras)
            except IndexError:
                new_ras = volume_background.get_center_point()
                x, y, background_matrix = volume_background.slice_volume(projection, new_ras)
                x1, y1, overlay_1_matrix = volume_overlay_1.slice_volume(projection, new_ras)
                x2, y2, overlay_2_matrix = volume_overlay_2.slice_volume(projection, new_ras)
                self.logger.info("The volume center point has been used for %s snapshot of %s, %s and %s.", projection,
                                 background_path, overlay_1_path, overlay_2_path)

            self.writer.write_3_matrices(x, y, background_matrix, x1, y1, overlay_1_matrix, x2, y2, overlay_2_matrix,
                                         self.generate_file_name(projection, snapshot_name))

    def overlap_surface_annotation(self, surface_path, annotation, snapshot_name=SNAPSHOT_NAME):
        annotation = self.annotation_io.read(annotation)
        surface = self.read_surface(surface_path, False)
        self.writer.write_surface_with_annotation(surface, annotation,
                                                  self.generate_file_name('surface_annotation', snapshot_name))

    def overlap_volume_surfaces(self, volume_background, surfaces_path, use_center_surface, use_cc_point,
                                snapshot_name=SNAPSHOT_NAME):
        volume = self.volume_io.read(volume_background)

        if use_cc_point:
            ras = self.generic_io.get_ras_coordinates(self.read_t1_affine_matrix())
        else:
            ras = volume.get_center_point()

        surfaces = [self.read_surface(os.path.expandvars(surface), use_center_surface) for surface in surfaces_path]

        for projection in PROJECTIONS:
            try:
                x, y, background_matrix = volume.slice_volume(projection, ras)
            except IndexError:
                ras = volume.get_center_point()
                x, y, background_matrix = volume.slice_volume(projection, ras)
                self.logger.info("The volume center point has been used for %s snapshot of %s and %s.", projection,
                                 volume_background, surfaces_path)

            clear_flag = True
            for surface_index, surface in enumerate(surfaces):
                surf_x_array, surf_y_array = surface.cut_by_plane(projection, ras)
                self.writer.write_matrix_and_surfaces(x, y, background_matrix, surf_x_array, surf_y_array,
                                                      surface_index, clear_flag)
                clear_flag = False
            self.writer.save_figure(self.generate_file_name(projection, snapshot_name))

    def show_aparc_aseg_with_new_values(self, aparc_aseg_volume_path, region_values_path, background_volume_path,
                                        use_cc_point, fs_to_conn_indices_mapping_path=FS_TO_CONN_INDICES_MAPPING_PATH,
                                        snapshot_name=SNAPSHOT_NAME):

        aparc_aseg_volume = self.volume_io.read(aparc_aseg_volume_path)

        fs_to_conn_indices_mapping = {}
        with open(fs_to_conn_indices_mapping_path) as fs_to_conn_indices_mapping_file:
            for line in fs_to_conn_indices_mapping_file:
                (key, label, val) = line.split()
                fs_to_conn_indices_mapping[float(key)] = float(val)

        conn_measure = numpy.zeros(len(fs_to_conn_indices_mapping))
        idx = 0
        with open(region_values_path) as region_values_file:
            for line in region_values_file:
                conn_measure[idx] = float(line.strip())
                idx += 1

        if use_cc_point:
            ras = self.generic_io.get_ras_coordinates(self.read_t1_affine_matrix())
        else:
            ras = aparc_aseg_volume.get_center_point()

        if background_volume_path != '':
            background_volume = self.volume_io.read(background_volume_path)

        for projection in PROJECTIONS:
            try:
                x_axis_coords, y_axis_coords, aparc_aseg_matrix = aparc_aseg_volume.slice_volume(projection, ras)
            except IndexError:
                new_ras = aparc_aseg_volume.get_center_point()
                x_axis_coords, y_axis_coords, aparc_aseg_matrix = aparc_aseg_volume.slice_volume(projection, new_ras)
                self.logger.info("The volume center point has been used for %s snapshot of %s.", projection,
                                 aparc_aseg_volume_path)

            for i in xrange(aparc_aseg_matrix.shape[0]):
                for j in xrange(aparc_aseg_matrix.shape[1]):
                    if aparc_aseg_matrix[i][j] > 0:
                        if fs_to_conn_indices_mapping.has_key(aparc_aseg_matrix[i][j]):
                            aparc_aseg_matrix[i][j] = conn_measure[
                                fs_to_conn_indices_mapping.get(aparc_aseg_matrix[i][j])]
                        else:
                            aparc_aseg_matrix[i][j] = -1

            if background_volume_path == '':
                self.writer.write_matrix(x_axis_coords, y_axis_coords, aparc_aseg_matrix,
                                         self.generate_file_name(projection, snapshot_name), 'hot')
            else:
                try:
                    bx_axis_coords, by_axis_coords, bvolume_matrix = background_volume.slice_volume(projection, ras)
                except IndexError:
                    new_ras = aparc_aseg_volume.get_center_point()
                    bx_axis_coords, by_axis_coords, bvolume_matrix = background_volume.slice_volume(projection, new_ras)
                    self.logger.info("The volume center point has been used for %s snapshot of %s and %s.", projection,
                                     aparc_aseg_volume_path, background_volume_path)

                self.writer.write_2_matrices(bx_axis_coords, by_axis_coords, bvolume_matrix, x_axis_coords,
                                             y_axis_coords,
                                             aparc_aseg_matrix,
                                             self.generate_file_name(projection, snapshot_name))
