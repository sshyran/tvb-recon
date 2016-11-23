# -*- coding: utf-8 -*-

import argparse
import os

from bnm.recon.logger import get_logger
from bnm.recon.qc.image.processor import ImageProcessor
from bnm.recon.qc.image.transformer import ImageTransformer
from bnm.recon.qc.model.constants import *

arg_1vol = "1vol"
arg_2vols = "2vols"
arg_3vols = "3vols"
arg_surf_annot = "surf_annot"
arg_vol_surf = "vol_surf"
arg_vol_white_pial = "vol_white_pial"


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate a BNM snapshot")
    subparsers = parser.add_subparsers(title='Sub Commands', dest='subcommand')

    parser.add_argument("--ras_transform", help="Applies RAS transformation on volumes", action="store_true")
    parser.add_argument("--center_surface", help="Centers surfaces using cras", action="store_true")

    subcommand_1_vol = subparsers.add_parser(arg_1vol, help='Display a single volume')
    subcommand_2_vols = subparsers.add_parser(arg_2vols, help='Display 2 volumes overlapped')
    subcommand_3_vols = subparsers.add_parser(arg_3vols, help='Display 3 volumes overlapped')
    subcommand_surf_annot = subparsers.add_parser(arg_surf_annot, help='Display a surface with annotations')
    subcommand_vol_surf = subparsers.add_parser(arg_vol_surf, help='Display a surface overlapped on a volume')
    subcommand_vol_2surf = subparsers.add_parser(arg_vol_white_pial,
                                                 help='Display white and pial surfaces over a volume')

    subcommand_1_vol.add_argument("volume")

    subcommand_2_vols.add_argument("background")
    subcommand_2_vols.add_argument("overlay")

    subcommand_3_vols.add_argument("background")
    subcommand_3_vols.add_argument("overlay1")
    subcommand_3_vols.add_argument("overlay2")

    subcommand_vol_surf.add_argument("background")
    subcommand_vol_surf.add_argument("surfaces_list", nargs='+')

    subcommand_surf_annot.add_argument("surface")
    subcommand_surf_annot.add_argument("annotation")

    subcommand_vol_2surf.add_argument("background")
    subcommand_vol_2surf.add_argument("-resampled_surface_name", default='')
    return parser.parse_args()


if __name__ == "__main__":
    logger = get_logger(__name__)

    args = parse_arguments()
    abs_path = os.path.abspath(os.path.dirname(__file__))

    snapshots_directory = os.environ[SNAPSHOTS_DIRECTORY_ENVIRON_VAR]
    if snapshots_directory is "":
        snapshots_directory = SNAPSHOTS_DIRECTORY
        logger.warning(
            "There is no value assigned to %s environment variable. The snapshots will be in %s directory by default.",
            SNAPSHOTS_DIRECTORY_ENVIRON_VAR, SNAPSHOTS_DIRECTORY)

    snapshot_count = int(os.environ.get(SNAPSHOT_NUMBER_ENVIRON_VAR, 0))

    imageTransformer = ImageTransformer(abs_path)
    imageTransformer.use_ras_transform = args.ras_transform
    imageTransformer.use_center_surface = args.center_surface

    imageProcessor = ImageProcessor(snapshots_directory=snapshots_directory, snapshot_count=snapshot_count)

    if args.subcommand == arg_1vol:
        volume_path = imageTransformer.transform_single_volume(os.path.expandvars(args.volume))
        imageProcessor.show_single_volume(os.path.expandvars(volume_path))

    elif args.subcommand == arg_2vols:
        background, overlay = imageTransformer.transform_2_volumes(os.path.expandvars(args.background),
                                                                   os.path.expandvars(args.overlay))
        imageProcessor.overlap_2_volumes(os.path.expandvars(background), os.path.expandvars(overlay))

    elif args.subcommand == arg_3vols:
        background, overlay1, overlay2 = imageTransformer.transform_3_volumes(os.path.expandvars(args.background),
                                                                              os.path.expandvars(args.overlay1),
                                                                              os.path.expandvars(args.overlay2))
        imageProcessor.overlap_3_volumes(os.path.expandvars(background), os.path.expandvars(overlay1),
                                         os.path.expandvars(overlay2))

    elif args.subcommand == arg_surf_annot:
        imageProcessor.overlap_surface_annotation(os.path.expandvars(args.surface), os.path.expandvars(args.annotation))

    elif args.subcommand == arg_vol_surf:
        background, surfaces_paths_list = imageTransformer.transform_volume_surfaces(
            os.path.expandvars(args.background), args.surfaces_list)
        imageProcessor.overlap_volume_surfaces(os.path.expandvars(background), surfaces_paths_list)

    elif args.subcommand == arg_vol_white_pial:
        surfaces_path = os.environ[SURFACES_DIRECTORY_ENVIRON_VAR]
        background, surfaces_paths_list = imageTransformer.transform_volume_white_pial(
            os.path.expandvars(args.background),
            os.path.expandvars(
                args.resampled_surface_name),
            os.path.expandvars(surfaces_path))
        imageProcessor.overlap_volume_surfaces(os.path.expandvars(background), os.path.expandvars(surfaces_paths_list))

    try:
        for created_file in imageTransformer.created_files:
            os.remove(created_file)
        os.rmdir(imageTransformer.converted_files_directory_path)
    except OSError:
        logger.error("Cannot delete files created at transform step")
