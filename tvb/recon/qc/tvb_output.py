# -*- coding: utf-8 -*-
import argparse
import os
import numpy
import shutil
from tvb.recon.algo.service.volume import VolumeService
from tvb.recon.io.factory import IOUtils
from tvb.recon.io.generic import GenericIO


def create_tvb_dataset(mri_direc: os.PathLike,
                       region_details_direc: os.PathLike,
                       weights_file: os.PathLike,
                       tracts_file: os.PathLike,
                       out_dir: os.PathLike,
                       bring_t1=False):
    weights_matrix = numpy.loadtxt(str(weights_file), dtype='i', delimiter=' ')
    weights_matrix += weights_matrix.T

    tracts_matrix = numpy.loadtxt(str(tracts_file), dtype='f', delimiter=' ')
    tracts_matrix += tracts_matrix.T

    is_cortical_rm = numpy.genfromtxt(os.path.join(region_details_direc, "cortical.txt"), usecols=[0], dtype='i')
    region_names = numpy.genfromtxt(os.path.join(region_details_direc, "centers.txt"), usecols=[0], dtype="str")
    region_centers = numpy.genfromtxt(os.path.join(region_details_direc, "centers.txt"), usecols=[1, 2, 3])
    region_areas = numpy.genfromtxt(os.path.join(region_details_direc, "areas.txt"), usecols=[0])
    region_orientations = numpy.genfromtxt(os.path.join(region_details_direc, "average_orientations.txt"),
                                           usecols=[0, 1, 2])
    rm_idx = numpy.genfromtxt(os.path.join(region_details_direc, "rm_to_aparc_aseg.txt"), usecols=[0, 1], dtype='i')
    rm_index_dict = dict(zip(rm_idx[:, 0], rm_idx[:, 1]))
    print(rm_index_dict)

    genericIO = GenericIO()
    genericIO.write_connectivity_zip(out_dir, weights_matrix, tracts_matrix, is_cortical_rm, region_names,
                                     region_centers, region_areas, region_orientations)

    aparc_aseg_file = os.path.join(mri_direc, "aparc+aseg.nii.gz")
    aparc_aseg_volume = IOUtils.read_volume(aparc_aseg_file)

    volume_service = VolumeService()
    aparc_aseg_cor_volume = volume_service.change_labels_of_aparc_aseg(aparc_aseg_volume, rm_index_dict,
                                                                       weights_matrix.shape[0])
    IOUtils.write_volume(os.path.join(out_dir, "aparc+aseg-cor.nii.gz"), aparc_aseg_cor_volume)

    if bring_t1:
        shutil.copy2(os.path.join(mri_direc, "T1.nii.gz"), out_dir)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Convert pipeline output to TVB format")
    parser.add_argument("-p", help="Call from Pegasus WMS", required=False, action="store_true")

    parser.add_argument("mri_dir")
    parser.add_argument("rm_details_dir")
    parser.add_argument("weights_file")
    parser.add_argument("tracts_file")
    parser.add_argument("output_dir")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    if args.p:
        create_tvb_dataset(
            args.mri_dir,
            args.rm_details_dir,
            args.weights_file,
            args.tracts_file,
            args.output_dir
        )

    else:
        create_tvb_dataset(
            args.mri_dir,
            args.rm_details_dir,
            args.weights_file,
            args.tracts_file,
            args.output_dir,
            True
        )
