# -*- coding: utf-8 -*-

import nibabel
import numpy
from nibabel.gifti import GiftiDataArray
from nibabel.gifti import GiftiImage
from nibabel.gifti import GiftiMetaData
from nibabel.gifti import giftiio
from nibabel.freesurfer.io import read_geometry, write_geometry
from bnm.recon.qc.model.surface import Surface
from bnm.recon.logger import get_logger


class ABCSurfaceParser(object):
    """
    This will define the behaviour needed for a surface parser.
    """

    def read(self, data_file):
        raise NotImplementedError()

    def write(self, surface_obj, file_path):
        raise NotImplementedError()

    def read_transformation_matrix_from_metadata(self, metadata):
        raise NotImplementedError()

    def write_transformation_matrix(self, metadata):
        raise NotImplementedError()


TRANSFORM_MATRIX_GIFTI_KEYS = [['VolGeomX_R', 'VolGeomY_R', 'VolGeomZ_R', 'VolGeomC_R'],
                               ['VolGeomX_A', 'VolGeomY_A', 'VolGeomZ_A', 'VolGeomC_A'],
                               ['VolGeomX_S', 'VolGeomY_S', 'VolGeomZ_S', 'VolGeomC_S']]


class GiftiSurfaceParser(ABCSurfaceParser):
    """
    This class reads content of GIFTI surface files
    """
    logger = get_logger(__name__)

    def read(self, data_file):
        gifti_image = giftiio.read(data_file)
        image_metadata = gifti_image.meta.metadata
        self.logger.info("From the file %s the extracted metadata is %s", data_file, image_metadata)

        data_arrays = gifti_image.darrays
        vertices = data_arrays[0].data
        triangles = data_arrays[1].data

        vol_geom_center_ras = [0, 0, 0]
        vertices_metadata = vertices.metadata
        self.logger.info("The metadata from vertices data array is %s", vertices_metadata)
        vertices_coord_system = vertices.coordsys
        self.logger.info(
            "The coordinate system transform matrix from vertices data array is %s", vertices_coord_system)
        triangles_metadata = triangles.metadata
        self.logger.info("The metadata from triangles data array is %s", triangles_metadata)

        # TODO review how and if we read this point
        vol_geom_center_ras[0] = float(vertices_metadata['VolGeomC_R'])
        vol_geom_center_ras[1] = float(vertices_metadata['VolGeomC_A'])
        vol_geom_center_ras[2] = float(vertices_metadata['VolGeomC_S'])

        return Surface(vertices, triangles, vol_geom_center_ras, image_metadata, vertices_metadata,
                       vertices_coord_system, triangles_metadata)

    def write(self, surface_obj, file_path):
        image_metadata = GiftiMetaData().from_dict(surface_obj.image_metadata)
        vertices_metadata = GiftiMetaData().from_dict(surface_obj.vertices_metadata)
        triangles_metadata = GiftiMetaData().from_dict(surface_obj.triangles_metadata)

        gifti_image = GiftiImage()
        gifti_image.set_metadata(image_metadata)

        data = GiftiDataArray(surface_obj.vertices, datatype='NIFTI_TYPE_FLOAT32', intent='NIFTI_INTENT_POINTSET')
        data.meta = vertices_metadata
        data.coordsys = surface_obj.vertices_coord_system
        gifti_image.add_gifti_data_array(data)

        data = GiftiDataArray(surface_obj.triangles, datatype='NIFTI_TYPE_INT32', intent='NIFTI_INTENT_TRIANGLE')
        data.meta = triangles_metadata
        data.coordsys = None
        gifti_image.add_gifti_data_array(data)

        nibabel.save(gifti_image, file_path)

    def read_transformation_matrix_from_metadata(self, image_metadata):
        matrix_from_metadata = [[0, 0, 0, 0] for _ in xrange(4)]

        for i in xrange(3):
            for j in xrange(4):
                matrix_from_metadata[i][j] = float(
                    image_metadata[TRANSFORM_MATRIX_GIFTI_KEYS[i][j]])

        matrix_from_metadata[3] = [0.0, 0.0, 0.0, 1.0]
        return matrix_from_metadata

    def write_transformation_matrix(self, image_metadata):
        # we can temporary write the identity matrix to gifti meta to avoid freeview rotations.

        identity_matrix = [[1.0, 0.0, 0.0, 0.0],
                           [0.0, 1.0, 0.0, 0.0],
                           [0.0, 0.0, 1.0, 0.0]]

        for i in xrange(3):
            for j in xrange(4):
                image_metadata[TRANSFORM_MATRIX_GIFTI_KEYS[i][j]] = str(identity_matrix[i][j])


TRANSFORM_MATRIX_FS_KEYS = ['xras', 'yras', 'zras', 'cras']


class FreesurferParser(ABCSurfaceParser):
    logger = get_logger(__name__)

    def read(self, surface_path):
        vertices, triangles, metadata = read_geometry(surface_path, read_metadata=True)
        self.logger.info("From the file %s the extracted metadata is %s", surface_path, metadata)

        # TODO make this better and this is not good if the cras of a centered surface can be read.
        if 'cras' in metadata:
            cras = metadata['cras']
            self.logger.info("The ras centering point for surface %s is %s", surface_path, cras)
        else:
            cras = [0, 0, 0]
            self.logger.warning("Could not read the ras centering point from surface %s header. The cras will be %s",
                                surface_path, cras)

        return Surface(vertices, triangles, cras, metadata)


    def write(self, surface, surface_path):
        write_geometry(filepath=surface_path, coords=surface.vertices, faces=surface.triangles,
                       volume_info=surface.get_main_metadata())


    def read_transformation_matrix_from_metadata(self, image_metadata):
        matrix_from_metadata = [[0, 0, 0, 0] for _ in xrange(4)] #or numpy.zeros((4,4))

        for i, fs_key in enumerate(TRANSFORM_MATRIX_FS_KEYS):
            for j in xrange(3):
                matrix_from_metadata[i][j] = image_metadata[fs_key][j]
        matrix_from_metadata[3][3] = 1
        matrix_from_metadata = numpy.transpose(matrix_from_metadata)
        return matrix_from_metadata

    def write_transformation_matrix(self, image_metadata):
        """
        We write the identity matrix to FS meta to avoid freeview rotations.
        :param image_metadata: meta to be corrected
        :return: image_metadata after change
        """
        identity_matrix = [[1.0, 0.0, 0.0],
                           [0.0, 1.0, 0.0],
                           [0.0, 0.0, 1.0],
                           [0.0, 0.0, 0.0]]
        for i, fs_key in enumerate(TRANSFORM_MATRIX_FS_KEYS):
            image_metadata[fs_key] = identity_matrix[i]
