# -*- coding: utf-8 -*-

import nibabel
from bnm.recon.logger import get_logger
from bnm.recon.qc.model.volume import Volume


class VolumeParser(object):
    """
    This class reads content of a NIFTI file and returns a Volume Object
    """

    logger = get_logger(__name__)

    def read(self, data_file):
        image = nibabel.load(data_file)
        header = image.header
        data = image.get_data()
        affine_matrix = image.affine
        self.logger.info("The affine matrix extracted from volume %s is %s" % (data_file, affine_matrix))

        return Volume(data, affine_matrix, header)
