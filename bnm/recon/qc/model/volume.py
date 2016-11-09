# -*- coding: utf-8 -*-

import numpy
import numpy.linalg as nlp
from nibabel.affines import apply_affine
from bnm.recon.qc.model.constants import *


class Volume(object):
    """
    Hold volume data, dimensions and affine matrix.

    Has a method that cuts a slice from the volume.
    """

    # 3D array
    data = None
    # array with the length of each data dimension
    dimensions = None
    # matrix containing voxel to ras transformation
    affine_matrix = None
    # TODO use this
    # transparency used when displaying a volume slice
    transparency = 1.0
    # color map used when displaying a volume slice
    color_map = None  # Reference to ColorMap obj, or empty

    def __init__(self, data, affine_matrix):
        self.data = data
        self.dimensions = data.shape
        self.affine_matrix = affine_matrix

    def slice_volume(self, projection=sagittal, ras=origin):
        """
        This determines slice colors and axes coordinates for the slice.
        :param projection: one of sagittal, axial or coronal
        :param ras: 3D point where to do the slicing
        :return: X, Y, 2D data matrix
        """

        affine_inverse = nlp.inv(self.affine_matrix)
        ijk_ras = numpy.round(apply_affine(affine_inverse, ras)).astype('i')

        slice_index_1, slice_index_2 = x_y_index[projection]

        slice_data = numpy.zeros((self.dimensions[slice_index_1], self.dimensions[slice_index_2]))
        x_axis_coords = numpy.zeros((self.dimensions[slice_index_1], self.dimensions[slice_index_2]))
        y_axis_coords = numpy.zeros((self.dimensions[slice_index_1], self.dimensions[slice_index_2]))

        for i in xrange(self.dimensions[slice_index_1]):
            for j in xrange(self.dimensions[slice_index_2]):
                ijk_ras[slice_index_1] = i
                ijk_ras[slice_index_2] = j

                ras_coordinates = apply_affine(self.affine_matrix, ijk_ras)
                x_axis_coords[i, j] = ras_coordinates[slice_index_1]
                y_axis_coords[i, j] = ras_coordinates[slice_index_2]

                color = self.data[ijk_ras[0], ijk_ras[1], ijk_ras[2]]
                if isinstance(color, (list, numpy.ndarray)):
                    color = color[0]
                slice_data[i][j] = color

        return x_axis_coords, y_axis_coords, slice_data
