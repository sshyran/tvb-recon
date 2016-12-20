# -*- coding: utf-8 -*-

import argparse
import os
import numpy
from bnm.recon.io.generic import GenericIO
from bnm.recon.logger import get_logger
from bnm.recon.model.constants import SNAPSHOTS_DIRECTORY_ENVIRON_VAR, SNAPSHOT_NUMBER_ENVIRON_VAR
from bnm.recon.qc.image.processor import ImageProcessor


def parse_arguments():
    parser = argparse.ArgumentParser(description="Transform a surface from its native space to another space")

    parser.add_argument("surface_path")
    parser.add_argument("output_path")
    parser.add_argument("-matrix_paths", nargs='+', default=[])
    parser.add_argument("-ss", help='Create snapshots of the transformed surface', action="store_true")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    surface_path = os.path.expandvars(args.surface_path)
    output_path = os.path.expandvars(args.output_path)

    logger = get_logger(__name__)

    image_processor = ImageProcessor(snapshots_directory=os.environ[SNAPSHOTS_DIRECTORY_ENVIRON_VAR],
                                     snapshot_count=int(os.environ.get(SNAPSHOT_NUMBER_ENVIRON_VAR, 0)))
    generic_io = GenericIO()

    logger.info("The surface transformation process has began")
    surface_parser = image_processor.factory_surface_io(surface_path)
    surface = surface_parser.read(surface_path, False)

    if len(args.matrix_paths) is not 0:
        transformation_matrices = []

        for transform_matrix_path in args.matrix_paths:
            transformation_matrices.append(numpy.array(generic_io.read_transformation_matrix(os.path.expandvars(transform_matrix_path))))

        for i in xrange(len(surface.vertices)):
            for j in xrange(len(transformation_matrices)):
                if len(transformation_matrices[j]) > 3:
                    vertex_coords = numpy.array(
                        [surface.vertices[i][0], surface.vertices[i][1], surface.vertices[i][2], 1])
                else:
                    vertex_coords = surface.vertices[i]

                new_vertex_coords = transformation_matrices[j].dot(vertex_coords)
                surface.vertices[i] = new_vertex_coords[:3]

    else:
        main_metadata = surface.get_main_metadata()
        transform_matrix = surface_parser.read_transformation_matrix_from_metadata(main_metadata)
        for i in xrange(len(surface.vertices)):
            vertex_coords = numpy.array([surface.vertices[i][0], surface.vertices[i][1], surface.vertices[i][2], 1])
            new_vertex_coords = vertex_coords.dot(transform_matrix)
            surface.vertices[i] = new_vertex_coords[:3]
        surface_parser.write_transformation_matrix(main_metadata)
        surface.set_main_metadata(main_metadata)

    surface_parser.write(surface, output_path)
    logger.info("The transformed surface has been written to file: %s" % output_path)

    if args.ss:
        image_processor.writer.write_surface(surface, image_processor.generate_file_name("transformed_surface"))
