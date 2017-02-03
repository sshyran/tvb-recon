# -*- coding: utf-8 -*-

import os
import time
from tvb.recon.logger import get_logger


def generate_uq_file_name(parent_folder, name_pattern):
    file_name = name_pattern.replace("*", str(time.clock()))
    return os.path.join(parent_folder, file_name)


def write_dax(adag, dax_file_name):
    # Write the DAX to stdout
    logger = get_logger(__name__)
    logger.info("Writing DAX into %s: %s" % (dax_file_name, adag))
    with open(dax_file_name, "w") as f:
        adag.writeXML(f)
