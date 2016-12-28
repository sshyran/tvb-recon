# -*- coding: utf-8 -*-

import numpy
from bnm.recon.algo.service.volume import VolumeService
from bnm.recon.model.volume import Volume


def _test_label_vol_from_tdi():
    service = VolumeService()
    data = numpy.array([[[0, 0, 1], [1, 2, 0]], [[2, 1, 3], [3, 1, 0]], [[0, 0, 1], [1, 2, 0]], [[2, 1, 3], [3, 1, 0]]])
    volume = Volume(data, [], None)
    labeled_volume = service._label_volume(volume, 0.5)
    assert labeled_volume.data.all() == numpy.array(
        [[[0, 0, 1], [2, 3, 0]], [[4, 5, 6], [7, 8, 0]], [[0, 0, 9], [10, 11, 0]], [[12, 13, 14], [15, 16, 0]]]).all()


def test_simple_label_config():
    service = VolumeService()
    data = numpy.array([[[0, 0, 1], [1, 2, 0]], [[2, 1, 3], [3, 1, 0]], [[0, 0, 1], [1, 2, 0]], [[2, 1, 3], [3, 1, 0]]])
    in_volume = Volume(data, [], None)
    out_volume = service._label_config(in_volume)
    assert out_volume.data.all() == numpy.array(
        [[[0, 0, 1], [1, 2, 0]], [[2, 1, 3], [3, 1, 0]], [[0, 0, 1], [1, 2, 0]], [[2, 1, 3], [3, 1, 0]]]).all()
