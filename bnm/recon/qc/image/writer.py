# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as pyplot
from mpl_toolkits.mplot3d import Axes3D
from bnm.recon.logger import get_logger
from bnm.recon.qc.model.constants import SNAPSHOT_EXTENSION


class ImageWriter(object):
    logger = get_logger(__name__)

    def __init__(self, snapshots_directory):
        self.snapshots_directory = snapshots_directory

        if not os.path.exists(self.snapshots_directory):
            os.mkdir(self.snapshots_directory)

    def get_path(self, result_name):
        return self.snapshots_directory + '/' + result_name + SNAPSHOT_EXTENSION

    def write_matrix(self, x, y, matrix, result_name):
        pyplot.pcolormesh(x, y, matrix, cmap="gray")
        pyplot.axis('off')
        pyplot.savefig(self.get_path(result_name), bbox_inches='tight', pad_inches=0.0)
        pyplot.clf()

    def write_2_matrices(self, x, y, matrix_background, x1, y1, matrix_overlap, result_name):
        pyplot.pcolormesh(x, y, matrix_background, cmap="gray")
        pyplot.pcolormesh(x1, y1, matrix_overlap, cmap="hot", alpha=0.3)
        pyplot.axis('off')
        pyplot.savefig(self.get_path(result_name), bbox_inches='tight', pad_inches=0.0)
        pyplot.clf()

    def write_3_matrices(self, x, y, matrix_background, x1, y1, matrix_overlap_1, x2, y2, matrix_overlap_2,
                         result_name):
        pyplot.pcolormesh(x, y, matrix_background, cmap="gray")
        pyplot.pcolormesh(x1, y1, matrix_overlap_1, cmap="hot", alpha=0.3)
        pyplot.pcolormesh(x2, y2, matrix_overlap_2, cmap="jet", alpha=0.5)
        pyplot.axis('off')
        pyplot.savefig(self.get_path(result_name), bbox_inches='tight', pad_inches=0.0)

    def write_surface(self, surface, result_name, positions=[(0, 0), (0, 90), (0, 180), (0, 270), (90, 0), (270, 0)]):
        #TODO show bigger surfaces in snapshots
        figs_folder = self.snapshots_directory
        self.logger.info("6 snapshots of the 3D surface will be generated in folder: %s" % figs_folder)

        x = surface.vertices[:, 0]
        y = surface.vertices[:, 1]
        z = surface.vertices[:, 2]

        fig = pyplot.figure()

        ax = Axes3D(fig)
        ax.set_xlim3d(-120, 60)
        ax.set_ylim3d(-120, 120)
        ax.set_zlim3d(-60, 120)
        ax.dist = 3

        ax.plot_trisurf(x, y, z, triangles=surface.triangles, color='g', linewidth=0.4)
        pyplot.axis('off')

        snapshot_index = 0
        for e, a in positions:
            ax.view_init(elev=e, azim=a)
            pyplot.savefig(self.get_path(result_name + str(snapshot_index)))
            snapshot_index += 1

        self.logger.info("The 6 snapshots were generated")

    def write_surface_with_annotation(self, surface, annot, result_name,
                                      positions=[(0, 0), (0, 90), (0, 180), (0, 270), (90, 0), (270, 0)]):
        x = surface.vertices[:, 0]
        y = surface.vertices[:, 1]
        z = surface.vertices[:, 2]

        fig = pyplot.figure()

        ax = Axes3D(fig)
        ax.set_xlim3d(-120, 60)
        ax.set_ylim3d(-120, 120)
        ax.set_zlim3d(-60, 120)
        ax.dist = 4

        face_colors = annot.compute_face_colors(surface.triangles)

        normals = surface.compute_normals()
        face_colors = ax._shade_colors(face_colors, normals)

        poly_line = ax.plot_trisurf(x, y, z, triangles=surface.triangles)
        poly_line.set_edgecolor(face_colors)
        poly_line.set_facecolor(face_colors)

        pyplot.axis('off')

        snapshot_index = 0
        for e, a in positions:
            ax.view_init(elev=e, azim=a)
            pyplot.savefig(self.get_path(result_name + str(snapshot_index)))
            snapshot_index += 1

    def save_figure(self, result_name):
        pyplot.axis('off')
        pyplot.savefig(self.get_path(result_name), bbox_inches='tight', pad_inches=0.0)

    def write_matrix_and_surfaces(self, x_axis_coords, y_axis_coords, matrix_background, surface_x_array, surface_y_array, clear_flag):
        if clear_flag:
            pyplot.clf()
        pyplot.pcolormesh(x_axis_coords, y_axis_coords, matrix_background, cmap="gray")
        for contour in xrange(len(surface_x_array)):
            pyplot.plot(surface_x_array[contour][:], surface_y_array[contour][:], 'y')
