# -*- coding: utf-8 -*-

SAGITTAL = 'sagittal'
CORONAL = 'coronal'
AXIAL = 'axial'

PROJECTIONS = [SAGITTAL, CORONAL, AXIAL]

PLANE_NORMALS = {
    SAGITTAL: (1, 0, 0),
    CORONAL: (0, 1, 0),
    AXIAL: (0, 0, 1)
}

X_Y_INDEX = {
    SAGITTAL: [1, 2],
    CORONAL: [0, 2],
    AXIAL: [0, 1]
}

ORIGIN = [0, 0, 0]
