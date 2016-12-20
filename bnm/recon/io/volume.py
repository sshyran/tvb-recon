# -*- coding: utf-8 -*-

import nibabel
from bnm.recon.logger import get_logger
from bnm.recon.qc.model.volume import Volume


class VolumeIO(object):
    """
    This class reads content of a NIFTI file and returns a Volume Object
    """

    logger = get_logger(__name__)

    def read(self, volume_path):
        image = nibabel.load(volume_path)
        header = image.header
        data = image.get_data()
        affine_matrix = image.affine
        self.logger.info("The affine matrix extracted from volume %s is %s" % (volume_path, affine_matrix))

        return Volume(data, affine_matrix, header)

    def write(self, out_volume_path, volume):
        image = nibabel.Nifti1Image(volume.data, volume.affine_matrix, volume.header)
        nibabel.save(image, out_volume_path)