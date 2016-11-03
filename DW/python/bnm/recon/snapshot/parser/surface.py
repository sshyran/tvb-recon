# -*- coding: utf-8 -*-
import nibabel
from nibabel.gifti import GiftiDataArray
from nibabel.gifti import GiftiImage
from nibabel.gifti import GiftiMetaData
from nibabel.gifti import giftiio
from nibabel.freesurfer.io import read_geometry
from nibabel.freesurfer.io import write_geometry
from bnm.recon.snapshot.model.surface import Surface
from bnm.recon.logger import get_logger


class SurfaceParser(object):
    """
    This class reads content of a NIFTI file and returns a Volume Object
    """

    def parse_gifti(self, data_file):
        gifti_image = giftiio.read(data_file)
        image_metadata = gifti_image.meta.metadata

        logger = get_logger(__name__)
        logger.info("From the file %s the extracted metadata is %s" % (data_file, image_metadata))

        data_arrays = gifti_image.darrays
        vertices = data_arrays[0].data
        triangles = data_arrays[1].data

        vol_geom_center_ras = [0, 0, 0]
        vertices_metadata = data_arrays[0].metadata
        logger.info("The metadata from vertices data array is %s" % vertices_metadata)
        vertices_coord_system = data_arrays[0].coordsys
        logger.info("The coordinate system transform matrix from vertices data array is %s" % vertices_coord_system)
        triangles_metadata = data_arrays[1].metadata
        logger.info("The metadata from triangles data array is %s" % triangles_metadata)

        # TODO we could read this values directly form the metadata
        vol_geom_center_ras[0] = float(data_arrays[0].metadata['VolGeomC_R'])
        vol_geom_center_ras[1] = float(data_arrays[0].metadata['VolGeomC_A'])
        vol_geom_center_ras[2] = float(data_arrays[0].metadata['VolGeomC_S'])

        return Surface(vertices, triangles, vol_geom_center_ras, image_metadata, vertices_metadata, vertices_coord_system, triangles_metadata)


    def parse_fs(self, surface_path):
        (vertices, triangles, metadata) = read_geometry(surface_path, read_metadata=True)
        cras = metadata['cras']

        logger = get_logger(__name__)
        logger.info("From the file %s the extracted metadata is %s" % (surface_path, metadata))

        return Surface(vertices, triangles, cras, metadata)


    def write_gifti(self, surface, surface_path):
        gifti_image = GiftiImage()

        data_array = [0 for _ in xrange(2)]
        data_array[0] = surface.vertices
        data_array[1] = surface.triangles

        image_metadata = GiftiMetaData().from_dict(surface.image_metadata)
        vertices_metadata = GiftiMetaData().from_dict(surface.vertices_metadata)
        triangles_metadata = GiftiMetaData().from_dict(surface.triangles_metadata)

        # TODO We currently write metadata of the old surface
        gifti_image.set_metadata(image_metadata)

        data = GiftiDataArray(data_array[0], datatype='NIFTI_TYPE_FLOAT32', intent='NIFTI_INTENT_POINTSET')
        data.meta = vertices_metadata
        data.coordsys = surface.vertices_coord_system
        gifti_image.add_gifti_data_array(data)
        data = GiftiDataArray(data_array[1], datatype='NIFTI_TYPE_INT32', intent='NIFTI_INTENT_TRIANGLE')
        data.meta = triangles_metadata
        data.coordsys = None
        gifti_image.add_gifti_data_array(data)

        nibabel.save(gifti_image, surface_path)

    def write_fs(self, surface, surface_path):
        write_geometry(filepath=surface_path, coords=surface.vertices, faces=surface.triangles, volume_info=surface.image_metadata)