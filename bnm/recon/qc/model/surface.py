# -*- coding: utf-8 -*-

import numpy
from trimesh import Trimesh, intersections
from bnm.recon.qc.model.constants import *


class Surface(object):
    """
    Hold a surface mesh (vertices and triangles).

    Has also few methods to read from this mesh (e.g. a contour cut).
    """

    # array of (x,y,z) tuples
    vertices = None
    # array of (v1, v2, v3) indices in vertices array
    triangles = None
    # color used when displaying this surface as a contour
    # TODO use this
    contour_color = 'r'
    # [x, y, z]
    center_ras = None

    def __init__(self, vertices, triangles, center_ras, generic_metadata, vertices_metadata=None,
                 vertices_coord_system=None, triangles_metadata=None):

        self.vertices = vertices
        self.triangles = triangles
        self.center_ras = center_ras

        self.generic_metadata = generic_metadata
        self.vertices_metadata = vertices_metadata
        self.triangles_metadata = triangles_metadata

        self.vertices_coord_system = vertices_coord_system


    def get_main_metadata(self):
        if self.vertices_metadata is not None:
            return self.vertices_metadata
        return self.generic_metadata

    def set_main_metadata(self, new_metadata):
        if self.vertices_metadata is not None:
            self.vertices_metadata = new_metadata
        else:
            self.generic_metadata = new_metadata


    def _get_plane_origin(self, ras):
        plane_origin = numpy.subtract(ras, self.center_ras)
        return list(plane_origin)


    def cut_by_plane(self, projection=sagittal, ras=origin):
        """
        :param projection:
        :param ras:
        :return: Y_array, X_array
        """
        mesh = Trimesh(self.vertices, self.triangles)
        contours = intersections.mesh_plane(mesh, plane_normals[projection], self._get_plane_origin(ras))
        x_array = [0 for _ in xrange(len(contours))]
        y_array = [0 for _ in xrange(len(contours))]

        for s in xrange(0, len(contours)):
            x_array[s] = contours[s][:, x_y_index[projection][0]]
            y_array[s] = contours[s][:, x_y_index[projection][1]]

        return x_array, y_array

    def compute_normals(self):
        """
        :return: array of triangle normal vectors
        """
        normals = [[0 for _ in xrange(0, 3)] for _ in xrange(0, len(self.triangles))]

        for i, tri in enumerate(self.triangles):
            u = self.vertices[tri[1]] - self.vertices[tri[0]]
            v = self.vertices[tri[2]] - self.vertices[tri[0]]
            normals[i][0] = u[1] * v[2] - u[2] * v[1]
            normals[i][1] = u[2] * v[0] - u[0] * v[2]
            normals[i][2] = u[0] * v[1] - u[1] * v[0]

        return normals
