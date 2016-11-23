# -*- coding: utf-8 -*-

import os


temporary_folder = 'data/temp'
created_test_files = []

def get_data_file(*args):
    file_path = os.path.join('data', *args)

    if not os.path.exists(file_path):
        raise Exception("File not found %s" % file_path)

    return file_path

def get_temporary_files_path(*args):
    if not os.path.exists(temporary_folder):
        os.makedirs(temporary_folder)

    file_path = os.path.join(temporary_folder, *args)
    created_test_files.append(file_path)

    return file_path

def remove_temporary_test_files():
    for created_test_file in created_test_files:
        os.remove(created_test_file)
    os.rmdir(temporary_folder)