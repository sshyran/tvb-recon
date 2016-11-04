# -*- coding: utf-8 -*-

from nibabel.affines import apply_affine
import numpy.linalg as nlp
import numpy

from ..model.constants import sagittal, coronal, axial


class Volume(object):
    dims = None
    ras_center_point = [128, 128, 128]
    data = None  # 3D matrix
    affine_matrix = None
    transparency = 1.0
    color_map = None  # Reference to ColorMap obj, or empty

    def __init__(self, nifti_data, dims, nifti_affine_matrix):
        self.data = nifti_data
        self.dims = dims
        self.affine_matrix = nifti_affine_matrix

    def get_axial(self, index):
        return self.data[:,index]

    def align(self, projection, ras):
        
        inv = nlp.inv(self.affine_matrix)
        ijk_ras = numpy.round(apply_affine(inv, ras)).astype('i')
        # TODO align diffusion images  which are smaller
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

        ras_new_center = (numpy.array(self.dims[:3])) / 2
        #ras_new_center = self.affine_matrix.dot(list(volume_center_voxel) + [1])
        #ras_new_center = map(int, ras_new_center)

        imin = - self.ras_center_point[ras_index_1] + ras_new_center[ras_index_1]
        imax = self.ras_center_point[ras_index_1] + ras_new_center[ras_index_1]
        jmin = - self.ras_center_point[ras_index_2] + ras_new_center[ras_index_2]
        jmax = self.ras_center_point[ras_index_2] + ras_new_center[ras_index_2]

        aligned_data = [[0 for _ in xrange(jmax-jmin)] for _ in xrange(imax-imin)]
        X=numpy.zeros((self.data.shape[ras_index_1],self.data.shape[ras_index_2]))
        Y=numpy.array(X)
        for i in range(imin, imax):
            for j in range(jmin, jmax):
                ras_current_point[ras_index_1] = i
                ras_current_point[ras_index_2] = j
                v = apply_affine(self.affine_matrix, ras_current_point)
                X[i,j]=v[ras_index_1]
                Y[i,j]=v[ras_index_2]
                #v = map(int, v)
                #if 0 <= v[0] < self.dims[0] and 0 <= v[1] < self.dims[1] and 0 <= v[2] < self.dims[2]:
                val = self.data[ras_current_point[0], ras_current_point[1], ras_current_point[2]]
                #if isinstance(val, (list, numpy.ndarray)):
                #val = val[0]
                            #else:
                            #val = 0
                aligned_data[i - imin][j - jmin] = val
        return (X,Y,aligned_data)

class ColorMap(object):
    # dict(VolumeValue : [r, g, b])
    colors = dict()
