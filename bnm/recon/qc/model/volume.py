# -*- coding: utf-8 -*-

from nibabel.affines import apply_affine
import numpy.linalg as nlp
import numpy

from bnm.recon.snapshot.model.constants import sagittal, coronal, axial


class Volume(object):
    data = None  # 3D matrix
    dims = None
    affine_matrix = None
    transparency = 1.0
    color_map = None  # Reference to ColorMap obj, or empty

    def __init__(self, nifti_data, dims, nifti_affine_matrix):
        self.data = nifti_data
        self.dims = dims
        self.affine_matrix = nifti_affine_matrix

    def align(self, projection, ras):

        inv = nlp.inv(self.affine_matrix)
        ijk_ras = numpy.round(apply_affine(inv, ras)).astype('i')
        ras_current_point = [0 for x in range(3)]

        if projection == sagittal:
            ras_current_point[0] = ijk_ras[0]
            ras_index_1 = 1
            ras_index_2 = 2
        elif projection == coronal:
            ras_current_point[1] = ijk_ras[1]
            ras_index_1 = 0
            ras_index_2 = 2
        elif projection == axial:
            ras_current_point[2] = ijk_ras[2]
            ras_index_1 = 0
            ras_index_2 = 1

        aligned_data = [[0 for _ in xrange(self.dims[ras_index_2])] for _ in xrange(self.dims[ras_index_1])]
        X = numpy.zeros((self.data.shape[ras_index_1], self.data.shape[ras_index_2]))
        Y = numpy.zeros((self.data.shape[ras_index_1], self.data.shape[ras_index_2]))

        for i in xrange(0, self.dims[ras_index_1]):
            for j in xrange(0, self.dims[ras_index_2]):
                ras_current_point[ras_index_1] = i
                ras_current_point[ras_index_2] = j

                v = apply_affine(self.affine_matrix, ras_current_point)
                X[i, j] = v[ras_index_1]
                Y[i, j] = v[ras_index_2]

                val = self.data[ras_current_point[0], ras_current_point[1], ras_current_point[2]]
                if isinstance(val, (list, numpy.ndarray)):
                    val = val[0]
                aligned_data[i][j] = val

        return X, Y, numpy.array(aligned_data)


class ColorMap(object):
    # dict(VolumeValue : [r, g, b])
    colors = dict()
