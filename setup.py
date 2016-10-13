# -*- coding: utf-8 -*-

"""
Prepare bnm Python package for setup.
"""

import shutil
from setuptools import setup, find_packages


setup(
    name='bnm',
    packages=find_packages(),
    version="1.0",
    author="BNM Team"
)

shutil.rmtree('bnm.egg-info', True)
