# -*- coding: utf-8 -*-

sagittal = 'sagittal'
coronal = 'coronal'
axial = 'axial'

projections = {sagittal, coronal, axial}

plane_normals = {
    sagittal: (1, 0, 0),
    coronal: (0, 1, 0),
    axial: (0, 0, 1)
}

x_y_index = {
    sagittal: [1, 2],
    coronal: [0, 2],
    axial: [0, 1]
}

origin = [0, 0, 0]
