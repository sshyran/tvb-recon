import numpy
import os
from tvb.recon.io.generic import GenericIO
from tvb.recon.io.volume import VolumeIO


def generate_schema_txt(ct_labeled_volume, schema_dir, schema_file):
    sensor_name = "sensor%s_"

    volumeIO = VolumeIO()
    vol = volumeIO.read(ct_labeled_volume)
    labels_nr = numpy.unique(vol.data).size

    schema_dict = dict()
    schema_dict.update({lbl: sensor_name % lbl for lbl in range(labels_nr)})

    genericIO = GenericIO()
    genericIO.write_dict_to_txt_file(schema_dict, os.path.join(schema_dir, schema_file))


def read_sensors_positions(sensors_file):
    sensors_positions = numpy.genfromtxt(sensors_file, usecols=[1, 2, 3])
    sensors_labels = numpy.genfromtxt(sensors_file, usecols=[0], dtype="str")

    return sensors_positions, sensors_labels
