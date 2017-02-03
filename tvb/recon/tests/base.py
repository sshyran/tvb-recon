# -*- coding: utf-8 -*-

import os
import glob


temporary_folder = 'data/temp'


def get_data_file(*args):
    file_path = os.path.join('data', *args)

    if not os.path.exists(file_path):
        raise Exception("File not found %s" % file_path)

    return file_path


def get_temporary_files_path(*args):
    if not os.path.exists(temporary_folder):
        os.makedirs(temporary_folder)

    file_path = os.path.join(temporary_folder, *args)
    return file_path


def remove_temporary_test_files():
    files = glob.glob('data/temp/*')
    for f in files:
        os.remove(f)
    os.rmdir(temporary_folder)
