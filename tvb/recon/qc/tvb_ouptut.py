# -*- coding: utf-8 -*-

import os

import sys

from tvb.recon.algo.service.annotation import AnnotationService
from tvb.recon.algo.service.surface import SurfaceService
from tvb.recon.io.factory import IOUtils


def create_tvb_dataset(subj_dir: os.PathLike,
                       cort_surf_direc: os.PathLike,
                       label_direc: os.PathLike,
                       out_surfaces_dir: os.PathLike = None):

    annotation_service = AnnotationService()
    cort_lut_file, subcort_lut_file = annotation_service.generate_cort_and_subcort_lut(subj_dir)

    annotation_io = IOUtils.annotation_io_factory("lh.aparc.annot")
    annot_cort_lh = annotation_io.read(os.path.join(label_direc, "lh.aparc.annot"))
    annot_cort_rh = annotation_io.read(os.path.join(label_direc, "rh.aparc.annot"))

    cort_region_mapping = annotation_service.generate_region_mapping(annot_cort_lh, annot_cort_rh, cort_lut_file, "aparc")

    annot_subcort_lh = annotation_io.read(os.path.join(label_direc, "lh.aseg.annot"))
    annot_subcort_rh = annotation_io.read(os.path.join(label_direc, "rh.aseg.annot"))

    subcort_region_mapping = annotation_service.generate_region_mapping(annot_subcort_lh, annot_subcort_rh,
                                                                        subcort_lut_file, "aseg")

    surface_service = SurfaceService()

    surface_io = IOUtils.surface_io_factory("lh.pial")
    surf_cort_lh = surface_io.read(os.path.join(cort_surf_direc, "lh.pial"), False)
    surf_cort_rh = surface_io.read(os.path.join(cort_surf_direc, "rh.pial"), False)

    full_cort_surface = surface_service.merge_surfaces([surf_cort_lh, surf_cort_rh])

    surf_subcort_lh = surface_io.read(os.path.join(cort_surf_direc, "lh.aseg"), False)
    surf_subcort_rh = surface_io.read(os.path.join(cort_surf_direc, "rh.aseg"), False)

    full_subcort_surface = surface_service.merge_surfaces([surf_subcort_lh, surf_subcort_rh])

    if out_surfaces_dir:
        with open(os.path.join(out_surfaces_dir, "region_mapping_cort.txt"), "w") as f:
            for rm_val in cort_region_mapping:
                f.write("%s\n" % rm_val)
        with open(os.path.join(out_surfaces_dir, "region_mapping_subcort.txt"), "w") as f:
            for rm_val in subcort_region_mapping:
                f.write("%s\n" % rm_val)
        io_utils = IOUtils()

        surf_subcort_filename = "surface_subcort.zip"
        surface_io = io_utils.surface_io_factory(surf_subcort_filename)
        surface_io.write(full_subcort_surface, os.path.join(out_surfaces_dir, surf_subcort_filename))

        surf_cort_filename = "surface_cort.zip"
        surface_io = io_utils.surface_io_factory(surf_cort_filename)
        surface_io.write(full_cort_surface, os.path.join(out_surfaces_dir, surf_cort_filename))


if __name__ == "__main__":
    subject_dir, out_surf = sys.argv[1:]

    create_tvb_dataset(
        subject_dir,
        os.path.join(subject_dir, "surf"),
        os.path.join(subject_dir, "label"),
        out_surf
    )
