# -*- coding: utf-8 -*-

from bnm.recon.qc.parser.generic import GenericParser
from bnm.tests.base import get_data_file


def test_read_cc_point():
    parser = GenericParser()
    file_path = get_data_file("cc_point.txt")
    cc_point = parser.read_cc_point(file_path, GenericParser.point_line_flag)
    assert cc_point == [0.0, 0.0, 0.0, 1]
