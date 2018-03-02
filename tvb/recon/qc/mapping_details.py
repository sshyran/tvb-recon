import argparse
import subprocess
import os
import numpy
from tvb.recon.algo.service.annotation import AnnotationService
from tvb.recon.algo.service.mapping_service import MappingService
from tvb.recon.algo.service.surface import SurfaceService
from tvb.recon.dax import AtlasSuffix
from tvb.recon.dax.mappings import AsegFilesPy3
from tvb.recon.io.factory import IOUtils
from tvb.recon.io.generic import GenericIO

genericIO = GenericIO()


def compute_region_details(atlas: AtlasSuffix, fs_color_lut: os.PathLike, t1: os.PathLike, lh_cort: os.PathLike,
                           rh_cort: os.PathLike,
                           lh_cort_annot: os.PathLike, rh_cort_annot: os.PathLike, lh_subcort: os.PathLike,
                           rh_subcort: os.PathLike, lh_subcort_annot: os.PathLike, rh_subcort_annot: os.PathLike):
    annot_cort_lh = IOUtils.read_annotation(lh_cort_annot)
    annot_cort_rh = IOUtils.read_annotation(rh_cort_annot)

    annot_subcort_lh = IOUtils.read_annotation(lh_subcort_annot)
    annot_subcort_rh = IOUtils.read_annotation(rh_subcort_annot)

    mapping = MappingService(annot_cort_lh, annot_cort_rh, annot_subcort_lh, annot_subcort_rh)
    mapping.generate_region_mapping_for_cort_annot(annot_cort_lh, annot_cort_rh)
    mapping.generate_region_mapping_for_subcort_annot(annot_subcort_lh, annot_subcort_rh)

    surface_service = SurfaceService()

    surf_cort_lh = IOUtils.read_surface(lh_cort, False)
    surf_cort_rh = IOUtils.read_surface(rh_cort, False)

    full_cort_surface = surface_service.merge_surfaces([surf_cort_lh, surf_cort_rh])

    surf_subcort_lh = IOUtils.read_surface(lh_subcort, False)
    surf_subcort_rh = IOUtils.read_surface(rh_subcort, False)

    full_subcort_surface = surface_service.merge_surfaces([surf_subcort_lh, surf_subcort_rh])

    genericIO.write_list_to_txt_file(mapping.cort_region_mapping, AsegFilesPy3.RM_CORT_TXT.format(atlas))
    genericIO.write_list_to_txt_file(mapping.subcort_region_mapping, AsegFilesPy3.RM_SUBCORT_TXT.format(atlas))

    vox2ras_file = "vox2ras.txt"
    subprocess.call(["mri_info", "--vox2ras", t1, "--o", vox2ras_file])

    surf_subcort_filename = "surface_subcort.zip"
    IOUtils.write_surface(surf_subcort_filename, full_subcort_surface)

    surf_cort_filename = "surface_cort.zip"
    IOUtils.write_surface(surf_cort_filename, full_cort_surface)

    os.remove(vox2ras_file)

    cort_subcort_full_surf = surface_service.merge_surfaces([full_cort_surface, full_subcort_surface])
    cort_subcort_full_region_mapping = mapping.cort_region_mapping + mapping.subcort_region_mapping

    dict_fs_custom = mapping.get_mapping_for_connectome_generation()
    genericIO.write_dict_to_txt_file(dict_fs_custom, AsegFilesPy3.FS_CUSTOM_TXT.format(atlas))

    region_areas = surface_service.compute_areas_for_regions(mapping.get_all_regions(), cort_subcort_full_surf,
                                                             cort_subcort_full_region_mapping)
    genericIO.write_list_to_txt_file(region_areas, AsegFilesPy3.AREAS_TXT.format(atlas))

    region_centers = surface_service.compute_centers_for_regions(mapping.get_all_regions(), cort_subcort_full_surf,
                                                                 cort_subcort_full_region_mapping)
    cort_subcort_lut = mapping.get_entire_lut()
    region_names = list(cort_subcort_lut.values())

    with open(AsegFilesPy3.CENTERS_TXT.format(atlas), "w") as f:
        for idx, (val_x, val_y, val_z) in enumerate(region_centers):
            f.write("%s %.2f %.2f %.2f\n" % (region_names[idx], val_x, val_y, val_z))

    region_orientations = surface_service.compute_orientations_for_regions(mapping.get_all_regions(),
                                                                           cort_subcort_full_surf,
                                                                           cort_subcort_full_region_mapping)

    lh_region_centers = surface_service.compute_centers_for_regions(mapping.get_lh_regions(), surf_cort_lh,
                                                                    mapping.lh_region_mapping)
    lh_region_orientations = surface_service.compute_orientations_for_regions(mapping.get_lh_regions(), surf_cort_lh,
                                                                              mapping.lh_region_mapping)
    with open(AsegFilesPy3.LH_DIPOLES_TXT.format(atlas), "w") as f:
        for idx, (val_x, val_y, val_z) in enumerate(lh_region_centers):
            f.write("%.2f %.2f %.2f %.2f %.2f %.2f\n" % (
                val_x, val_y, val_z, lh_region_orientations[idx][0], lh_region_orientations[idx][1],
                lh_region_orientations[idx][2]))

    rh_region_centers = surface_service.compute_centers_for_regions(mapping.get_rh_regions(), surf_cort_rh,
                                                                    mapping.rh_region_mapping)
    rh_region_orientations = surface_service.compute_orientations_for_regions(mapping.get_rh_regions(), surf_cort_rh,
                                                                              mapping.rh_region_mapping)
    with open(AsegFilesPy3.RH_DIPOLES_TXT.format(atlas), "w") as f:
        for idx, (val_x, val_y, val_z) in enumerate(rh_region_centers):
            f.write("%.2f %.2f %.2f %.2f %.2f %.2f\n" % (
                val_x, val_y, val_z, rh_region_orientations[idx][0], rh_region_orientations[idx][1],
                rh_region_orientations[idx][2]))

    numpy.savetxt(AsegFilesPy3.ORIENTATIONS_TXT.format(atlas), region_orientations, fmt='%.2f %.2f %.2f')

    annotation_service = AnnotationService()
    lut_dict, _, _ = annotation_service.read_lut(fs_color_lut, "name")
    rm_index_dict = mapping.get_mapping_for_aparc_aseg(lut_dict)
    genericIO.write_dict_to_txt_file(rm_index_dict, AsegFilesPy3.RM_TO_APARC_ASEG_TXT.format(atlas))

    genericIO.write_list_to_txt_file(mapping.is_cortical_region_mapping(), AsegFilesPy3.CORTICAL_TXT.format(atlas))


def parse_arguments():
    parser = argparse.ArgumentParser(description="Convert pipeline output to TVB format")
    parser.add_argument("-p", help="Call from Pegasus WMS", required=False, action="store_true")

    parser.add_argument("atlas_suffix")
    parser.add_argument("fs_color_lut")
    parser.add_argument("t1")
    parser.add_argument("lh_cort")
    parser.add_argument("rh_cort")
    parser.add_argument("lh_cort_annot")
    parser.add_argument("rh_cort_annot")
    parser.add_argument("lh_subcort")
    parser.add_argument("rh_subcort")
    parser.add_argument("lh_subcort_annot")
    parser.add_argument("rh_subcort_annot")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    if args.p:
        compute_region_details(args.atlas_suffix, args.fs_color_lut, args.t1, args.lh_cort, args.rh_cort,
                               args.lh_cort_annot, args.rh_cort_annot, args.lh_subcort, args.rh_subcort,
                               args.lh_subcort_annot, args.rh_subcort_annot)
