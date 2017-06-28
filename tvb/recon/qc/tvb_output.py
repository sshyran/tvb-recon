# -*- coding: utf-8 -*-

import os
import subprocess
import sys

import numpy
import shutil
from tvb.recon.algo.service.annotation import AnnotationService
from tvb.recon.algo.service.mapping_service import MappingService
from tvb.recon.algo.service.surface import SurfaceService
from tvb.recon.algo.service.volume import VolumeService
from tvb.recon.io.factory import IOUtils
from tvb.recon.io.generic import GenericIO


def create_tvb_dataset(cort_surf_direc: os.PathLike,
                       label_direc: os.PathLike,
                       mri_direc: os.PathLike,
                       weights_file: os.PathLike,
                       tracts_file: os.PathLike,
                       fs_color_lut: os.PathLike,
                       out_dir: os.PathLike):
    annotation_io = IOUtils.annotation_io_factory("lh.aparc.annot")
    annot_cort_lh = annotation_io.read(os.path.join(label_direc, "lh.aparc.annot"))
    annot_cort_rh = annotation_io.read(os.path.join(label_direc, "rh.aparc.annot"))

    annot_subcort_lh = annotation_io.read(os.path.join(label_direc, "lh.aseg.annot"))
    annot_subcort_rh = annotation_io.read(os.path.join(label_direc, "rh.aseg.annot"))

    mapping = MappingService(annot_cort_lh, annot_cort_rh, annot_subcort_lh, annot_subcort_rh)
    mapping.generate_region_mapping_for_cort_annot(annot_cort_lh, annot_cort_rh)
    mapping.generate_region_mapping_for_subcort_annot(annot_subcort_lh, annot_subcort_rh)

    surface_service = SurfaceService()

    surface_io = IOUtils.surface_io_factory("lh-centered.pial")
    surf_cort_lh = surface_io.read(os.path.join(cort_surf_direc, "lh-centered.pial"), False)
    surf_cort_rh = surface_io.read(os.path.join(cort_surf_direc, "rh-centered.pial"), False)

    full_cort_surface = surface_service.merge_surfaces([surf_cort_lh, surf_cort_rh])

    surf_subcort_lh = surface_io.read(os.path.join(cort_surf_direc, "lh-centered.aseg"), False)
    surf_subcort_rh = surface_io.read(os.path.join(cort_surf_direc, "rh-centered.aseg"), False)

    full_subcort_surface = surface_service.merge_surfaces([surf_subcort_lh, surf_subcort_rh])

    with open(str(os.path.join(out_dir, "region_mapping_cort.txt")), "w") as f:
        for rm_val in mapping.cort_region_mapping:
            f.write("%s\n" % rm_val)
    with open(str(os.path.join(out_dir, "region_mapping_subcort.txt")), "w") as f:
        for rm_val in mapping.subcort_region_mapping:
            f.write("%s\n" % rm_val)

    vox2ras_file = os.path.join(out_dir, "vox2ras.txt")
    subprocess.call(["mri_info", "--vox2ras", os.path.join(mri_direc, "T1.nii.gz"), "--o", vox2ras_file])

    io_utils = IOUtils()
    surf_subcort_filename = "surface_subcort.zip"
    surface_io = io_utils.surface_io_factory(surf_subcort_filename)
    surface_io.write(full_subcort_surface, os.path.join(out_dir, surf_subcort_filename))

    surf_cort_filename = "surface_cort.zip"
    surface_io = io_utils.surface_io_factory(surf_cort_filename)
    surface_io.write(full_cort_surface, os.path.join(out_dir, surf_cort_filename))

    os.remove(vox2ras_file)

    cort_subcort_full_surf = surface_service.merge_surfaces([full_cort_surface, full_subcort_surface])
    cort_subcort_full_region_mapping = mapping.cort_region_mapping + mapping.subcort_region_mapping

    region_areas = surface_service.compute_areas_for_regions(mapping.get_all_regions(),
                                                             cort_subcort_full_surf,
                                                             cort_subcort_full_region_mapping)

    region_centers = surface_service.compute_centers_for_regions(mapping.get_all_regions(), cort_subcort_full_surf,
                                                                 cort_subcort_full_region_mapping)

    cort_subcort_lut = mapping.get_entire_lut()
    region_names = list(cort_subcort_lut.values())

    region_orientations = surface_service.compute_orientations_for_regions(mapping.get_all_regions(),
                                                                           cort_subcort_full_surf,
                                                                           cort_subcort_full_region_mapping)

    weights_matrix = numpy.loadtxt(str(weights_file), dtype='i', delimiter=' ')
    weights_matrix += weights_matrix.T

    tracts_matrix = numpy.loadtxt(str(tracts_file), dtype='f', delimiter=' ')
    tracts_matrix += tracts_matrix.T

    genericIO = GenericIO()
    genericIO.write_connectivity_zip(out_dir, weights_matrix, tracts_matrix,
                                     mapping.is_cortical_region_mapping(), region_names, region_centers,
                                     region_areas, region_orientations)

    annotation_service = AnnotationService()
    lut_dict, _, _ = annotation_service.read_lut(fs_color_lut, "name")
    rm_index_dict = mapping.get_mapping_for_aparc_aseg(lut_dict)

    aparc_aseg_file = os.path.join(mri_direc, "aparc+aseg.nii.gz")
    volume_io = IOUtils.volume_io_factory(aparc_aseg_file)
    aparc_aseg_volume = volume_io.read(aparc_aseg_file)

    volume_service = VolumeService()
    aparc_aseg_cor_volume = volume_service.change_labels_of_aparc_aseg(aparc_aseg_volume, rm_index_dict,
                                                                       weights_matrix.shape[0])
    volume_io.write(os.path.join(out_dir, "aparc+aseg-cor.nii.gz"), aparc_aseg_cor_volume)

    shutil.copy2(os.path.join(mri_direc, "T1.nii.gz"), out_dir)


if __name__ == "__main__":
    subject_dir, weights_file, tracts_file, fs_color_lut, out_dir = sys.argv[1:]

    create_tvb_dataset(
        os.path.join(subject_dir, "surf"),
        os.path.join(subject_dir, "label"),
        os.path.join(subject_dir, "mri"),
        weights_file,
        tracts_file,
        fs_color_lut,
        out_dir
    )
