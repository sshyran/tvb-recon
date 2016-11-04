# -*- coding: utf-8 -*-

import numpy
import trimesh
from trimesh import intersections
from .constants import sagittal, coronal, axial


class Surface(object):

    vertices = None
    triangles = None
    center_value = 128
    vol_geom_center_ras = None
    gii_reference = None
    contour_color = None
    show_all = False

    plane_normals = {
        sagittal : (1,0,0),
        coronal : (0,1,0),
        axial : (0,0,1)
    }

    x_y_index = {
        sagittal : [1, 2],
        coronal : [0, 2],
        axial : [0, 1]
    }

    gifti_transform_matrix_key_words = [['VolGeomX_R', 'VolGeomY_R', 'VolGeomZ_R', 'VolGeomC_R'],
                                        ['VolGeomX_A', 'VolGeomY_A', 'VolGeomZ_A', 'VolGeomC_A'],
                                        ['VolGeomX_S', 'VolGeomY_S', 'VolGeomZ_S', 'VolGeomC_S']]

    fs_transform_matrix_key_words = ['xras', 'yras', 'zras', 'cras']

    def read_matrix_from_metadata(self, is_gifti):
        matrix_from_metadata = [[0 for _ in range(4)] for _ in range(4)]

        if is_gifti:
            for i in range(3):
                for j in range(4):
                    matrix_from_metadata[i][j] = float(self.vertices_metadata[self.gifti_transform_matrix_key_words[i][j]])
            matrix_from_metadata[3] = [0.0, 0.0, 0.0, 1.0]
        else:
            for i in range(len(self.fs_transform_matrix_key_words)):
                for j in range(3):
                    matrix_from_metadata[i][j] = self.image_metadata[self.fs_transform_matrix_key_words[i]][j]
            matrix_from_metadata[3][3] = 1
            matrix_from_metadata = numpy.transpose(matrix_from_metadata)
        return matrix_from_metadata


    def write_matrix_from_metadata(self, is_gifti):
        # we can temporary write the identity matrix to gifti meta to avoid freeview rotations.
        identity_matrix = [ [1.0, 0.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0, 0.0],
                            [0.0, 0.0, 1.0, 0.0]]

        if is_gifti:
            for i in range(3):
                for j in range(4):
                    self.vertices_metadata[self.gifti_transform_matrix_key_words[i][j]] = str(identity_matrix[i][j])
        else:
            identity_matrix = numpy.transpose(identity_matrix)
            for i in range(len(self.fs_transform_matrix_key_words)):
                self.image_metadata[self.fs_transform_matrix_key_words[i]] = identity_matrix[i]

    # This method is not currently used.
    # Applies the rotation_matrix on every surface vertex.
    # Rotates surface contour and it is displayed similar to freeview.
    def apply_rotation_matrix(self, contour):
        rotation_matrix = self.read_matrix_from_metadata()
        new_verts = [[0 for _ in range(3)] for _ in range(len(contour))]

        for i in range(0, len(contour)):
            new_verts[i] = rotation_matrix.dot(contour[i])

        return numpy.array(new_verts)


    def __init__(self, vertices, triangles, vol_geom_center_ras, image_metadata, vertices_metadata=None, vertices_coord_system=None, triangles_metadata=None):
        self.vertices = vertices
        self.triangles = triangles
        self.vol_geom_center_ras = vol_geom_center_ras
        self.image_metadata = image_metadata
        self.vertices_metadata = vertices_metadata
        self.vertices_coord_system = vertices_coord_system
        self.triangles_metadata = triangles_metadata


    def get_plane_origin(self, ras):
        plane_origin = numpy.subtract(ras, self.vol_geom_center_ras)
        return list(plane_origin)


    def get_x_y_array(self, projection, ras):
        mesh = trimesh.Trimesh(self.vertices, self.triangles)
        contours = intersections.mesh_plane(mesh, self.plane_normals[projection], self.get_plane_origin(ras))
        x_array = [0 for _ in range(len(contours))]
        y_array = [0 for _ in range(len(contours))]

        for s in range(0, len(contours)):
            x_array[s] = contours[s][:, self.x_y_index[projection][0]]
            y_array[s] = contours[s][:, self.x_y_index[projection][1]]

        return x_array, y_array


    def compute_normals(self):
        normals = [[0 for _ in range(0, 3)] for _ in range(0, len(self.triangles))]

        for i, tri in enumerate(self.triangles):
            u = self.vertices[tri[1]] - self.vertices[tri[0]]
            v = self.vertices[tri[2]] - self.vertices[tri[0]]
            normals[i][0] = u[1] * v[2] - u[2] * v[1]
            normals[i][1] = u[2] * v[0] - u[0] * v[2]
            normals[i][2] = u[0] * v[1] - u[1] * v[0]

        return normals