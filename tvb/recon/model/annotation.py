# -*- coding: utf-8 -*-

from typing import Union, Optional
import numpy


class Annotation(object):
    """
    Hold annotation information as region mapping, color mapping and names of regions.

    Has a method to compute face colors using vertices_color_mapping .
    """

    def __init__(self, region_mapping: Union[numpy.ndarray, list], regions_color_table: numpy.ndarray,
                 region_names: list):
        if len(region_mapping) == 0:
            self.region_mapping = numpy.empty((0,), dtype='i')
        else:
            # ndarray of region_names indices
            self.region_mapping = numpy.array(region_mapping)
        if len(regions_color_table) == 0:
            self.regions_color_table = numpy.empty((0, 5), dtype='i')
        else:
            # ndarray matrix that contains a rgba color array for each region
            self.regions_color_table = numpy.array(regions_color_table)
        self.region_names = region_names  # list of region names

    def set_region_mapping(self, new_region_mapping: numpy.ndarray):
        self.region_mapping = new_region_mapping

    def add_region_names_and_colors(self, new_region_names: list, new_region_colors: numpy.ndarray):
        # TODO: check if the following line is correct! I changed it from .append() to +=
        self.region_names += new_region_names
        self.regions_color_table = numpy.concatenate(
            (self.regions_color_table, new_region_colors), axis=0).astype('i')

    def add_region_mapping(self, new_region_mapping: Union[numpy.ndarray, list]):
        self.region_mapping = numpy.r_[self.region_mapping, new_region_mapping]

    # def stack_region_mapping(self):
    #     self.region_mapping = numpy.hstack(self.region_mapping)

    def get_region_mapping_by_indices(self, indices: Union[numpy.ndarray, list]):
        return self.region_mapping[indices]

    def compute_face_colors(self, triangles: numpy.ndarray) -> list:
        face_colors = []
        converted_colors = self.regions_color_table[:, :4] / 255.0

        for triangle in triangles:
            vertex1_color = converted_colors[
                self.region_mapping[triangle[0]]] * 0.33
            vertex2_color = converted_colors[
                self.region_mapping[triangle[1]]] * 0.33
            vertex3_color = converted_colors[
                self.region_mapping[triangle[2]]] * 0.33
            face_colors.append(vertex1_color + vertex2_color + vertex3_color)

        return face_colors
