# -*- coding: utf-8 -*-

import os


def get_data_file(*args):
    file_path = os.path.join('data', *args)

    if not os.path.exists(file_path):
        raise Exception("File not found %s" % file_path)

    return file_path
