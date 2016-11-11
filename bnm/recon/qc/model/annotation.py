# -*- coding: utf-8 -*-

import numpy


class Annotation(object):
    """
    Hold annotation information as region mapping, color mapping and names of regions.

    Has a method to compute face colors using vertices_color_mapping .
    """

    # array of region_names indices
    region_mapping = None
    # matrix that contains a rgba color array for each vertex
    regions_color_table = None
    # list of region names
    region_names = None

    def __init__(self, region_mapping, regions_color_table, region_names):
        self.region_mapping = region_mapping
        self.regions_color_table = regions_color_table
        self.region_names = region_names

    def _convert_rgba(self):
        """
        This converts the rgba values from [0, 255] to [0, 1] interval.
        :return: the converted matrix
        """
        converted_color_table = [[0 for _ in xrange(4)] for _ in xrange(len(self.regions_color_table))]

        for i in xrange(len(self.regions_color_table)):
            for j in xrange(4):
                converted_color_table[i][j] = round((self.regions_color_table[i][j] / 255.0), 2)

        return converted_color_table

    def compute_face_colors(self, triangles):
        face_colors = [0 for _ in xrange(len(triangles))]
        converted_colors = self._convert_rgba()

        for i, triangle in enumerate(triangles):
            vertex1_color = numpy.dot(converted_colors[self.region_mapping[triangle[0]]], 0.33)
            vertex2_color = numpy.dot(converted_colors[self.region_mapping[triangle[1]]], 0.33)
            vertex3_color = numpy.dot(converted_colors[self.region_mapping[triangle[2]]], 0.33)

            intermediate_color = numpy.add(vertex1_color, vertex2_color)
            face_colors[i] = numpy.add(intermediate_color, vertex3_color)

        return face_colors
