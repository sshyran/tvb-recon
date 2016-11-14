# -*- coding: utf-8 -*-

import argparse
import os
from bnm.recon.qc.image.processor import ImageProcessor
from bnm.recon.qc.image.transformer import ImageTransformer

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
    args = parse_arguments()
    abs_path = os.path.abspath(os.path.dirname(__file__))
    imageTransformer = ImageTransformer(abs_path)

    imageTransformer.use_ras_transform = args.ras_transform
    imageTransformer.use_center_surface = args.center_surface

    imageProcessor = ImageProcessor(snapshots_directory=os.environ['FIGS'],
                                    snapshot_count=int(os.environ.get('SNAPSHOT_NUMBER', 0)))

    if args.subcommand == arg_1vol:
        volume_path = imageTransformer.transform_single_volume(os.path.expandvars(args.volume))
        imageProcessor.show_single_volume(os.path.expandvars(volume_path))

    elif args.subcommand == arg_2vols:
        bkg, ovr = imageTransformer.transform_2_volumes(os.path.expandvars(args.background),
                                                        os.path.expandvars(args.overlay))
        imageProcessor.overlap_2_volumes(os.path.expandvars(bkg), os.path.expandvars(ovr))

    elif args.subcommand == arg_3vols:
        bkg, ovr1, ovr2 = imageTransformer.transform_3_volumes(os.path.expandvars(args.background),
                                                               os.path.expandvars(args.overlay1),
                                                               os.path.expandvars(args.overlay2))
        imageProcessor.overlap_3_volumes(os.path.expandvars(bkg), os.path.expandvars(ovr1), os.path.expandvars(ovr2))

    elif args.subcommand == arg_surf_annot:
        imageProcessor.overlap_surface_annotation(os.path.expandvars(args.surface), os.path.expandvars(args.annotation))

    elif args.subcommand == arg_vol_surf:
        bkg, surf = imageTransformer.transform_volume_surfaces(os.path.expandvars(args.background),
                                                               os.path.expandvars(args.surfaces_list))
        imageProcessor.overlap_volume_surface(os.path.expandvars(bkg), os.path.expandvars(surf))

    elif args.subcommand == arg_vol_white_pial:
        # TODO for transformations, have a separate folder with the white+pial surfaces
        imageProcessor.overlap_volume_surfaces(os.path.expandvars(args.background),
                                               os.path.expandvars(args.resampled_surface_name),
                                               surfaces_path=os.path.expandvars(os.environ['SURF']))

    try:
        for i in imageTransformer.created_files:
            os.remove(i)
        os.rmdir(imageTransformer.converted_files_directory_path)
    except OSError:
        print "Cannot delete files"
