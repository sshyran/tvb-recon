# -*- coding: utf-8 -*-

from bnm.recon.qc.parser.surface import FreesurferParser
from bnm.tests.base import get_data_file


def test_parse_fs_surface():
    parser = FreesurferParser()
    file_path = get_data_file("freesurfer_fsaverage", "surf", "lh.pial")
    surf = parser.read(file_path)
    assert len(surf.triangles) == 327680
