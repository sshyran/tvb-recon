# -*- coding: utf-8 -*-

import os
import numpy
import scipy
import gdist
from nibabel.freesurfer.io import write_geometry
from sklearn.metrics.pairwise import paired_distances
from scipy.sparse import csr_matrix
from bnm.recon.algo.service.annotation import AnnotationService
from bnm.recon.qc.parser.volume import VolumeIO
from bnm.recon.qc.parser.surface import FreesurferIO


class SurfaceService(object):
    def __init__(self):
        self.annotation_service = AnnotationService()
        self.surface_parser = FreesurferIO()

    def tri_area(self, tri):
        i, j, k = numpy.transpose(tri, (1, 0, 2))
        ij = j - i
        ik = k - i
        return numpy.sqrt(numpy.sum(numpy.cross(ij, ik) ** 2, axis=1)) / 2.0

    def compute_gdist_mat(self, surf_name='pial', max_distance=40.0):
        max_distance = float(max_distance)  # in case passed from sys.argv
        for h in 'rl':
            subjects_dir = os.environ['SUBJECTS_DIR']
            subject = os.environ['SUBJECT']
            surf_path = '%s/%s/surf/%sh.%s' % (subjects_dir, subject, h, surf_name)
            surface = self.surface_parser.read(surf_path, False)
            mat_path = '%s/%s/surf/%sh.%s.gdist.mat' % (subjects_dir, subject, h, surf_name)
            mat = gdist.local_gdist_matrix(surface.vertices, surface.triangles.astype('<i4'), max_distance=max_distance)
            scipy.io.savemat(mat_path, {'gdist': mat})

    def extract_subsurf(self, verts, faces, verts_mask):
        # verts and faces to keep
        verts_out = verts[verts_mask, :]
        face_mask = numpy.c_[verts_mask[faces[:, 0]], verts_mask[faces[:, 1]], verts_mask[faces[:, 2]]].all(axis=1)
        faces_out = faces[face_mask]

        # ...but the old vertices' indexes of faces have to be transformed to the new verts_out_inds:
        verts_out_inds, = numpy.where(verts_mask)
        for face_idx in range(faces_out.shape[0]):
            for vertex_idx in range(3):
                faces_out[face_idx, vertex_idx], = numpy.where(faces_out[face_idx, vertex_idx] == verts_out_inds)

        return verts_out, faces_out

    # Concatenate surfaces of specific labels to create a single annotated surface
    def aseg_surf_conc_annot(self, surf_path, out_surf_path, annot_path, label_indices,
                             lut_path=os.path.join(os.environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt')):

        label_names, color_table = self.annotation_service.lut_to_annot_names_ctab(lut_path=lut_path,
                                                                                   labels=label_indices)
        label_indices = numpy.array(label_indices.split()).astype('i')

        out_verts = []
        out_faces = []
        out_labels = []
        out_names = []
        out_color_table = []
        label_index = -1
        verts_number = 0
        volume_info = None

        for lbl in label_indices:
            l = int(lbl)
            this_surf_path = surf_path + "-%06d" % l

            if os.path.exists(this_surf_path):
                ind_l, = numpy.where(label_indices == lbl)
                out_names.append(label_names[ind_l])
                out_color_table.append(color_table[ind_l, :])
                label_index += 1
                surface = self.surface_parser.read(this_surf_path, False)
                faces = surface.triangles
                volume_info = surface.get_main_metadata()
                faces = faces + verts_number  # Update vertices indexes
                verts_number += surface.vertices.shape[0]
                out_verts.append(surface.vertices)
                out_faces.append(faces)
                out_labels.append(label_index * numpy.ones((surface.vertices.shape[0],), dtype='int64'))

        out_color_table = numpy.squeeze(numpy.array(out_color_table).astype('i'))
        out_verts = numpy.vstack(out_verts)
        out_faces = numpy.vstack(out_faces)
        out_labels = numpy.hstack(out_labels)

        write_geometry(out_surf_path, out_verts, out_faces, volume_info=volume_info)
        self.annotation_service.annotationIO.write(annot_path, out_labels, out_color_table, out_names)

    # It returns a sparse matrix of the connectivity among the vertices of a surface
    # mode: "sparse" (default) or "2D"
    def vertex_connectivity(self, v, f, mode="sparse", metric=None):
        # Get all pairs of vertex indexes that appear in each face
        f = numpy.r_[f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]
        # Remove repetitions
        f = numpy.vstack(set(map(tuple, f)))
        # Mark all existing pairs to 1
        n_v = v.shape[0]
        n_f = f.shape[0]
        if metric is None:
            con = csr_matrix((numpy.ones((n_f,)), (f[:, 0], f[:, 1])), shape=(n_v, n_v))
            if mode != "sparse":
                # Create non-sparse matrix
                con = con.todense()
        else:
            d = paired_distances(v[f[:, 0]], v[f[:, 1]], metric)
            if mode == "sparse":
                # Create sparse matrix
                con = csr_matrix((d, (f[:, 0], f[:, 1])), shape=(n_v, n_v))
        return con

    def __prepare_grid(self, vertex_neighbourhood):
        # Prepare grid if needed for possible use:
        if vertex_neighbourhood > 0:
            grid = numpy.meshgrid(range(-vertex_neighbourhood, vertex_neighbourhood + 1, 1),
                                  range(-vertex_neighbourhood, vertex_neighbourhood + 1, 1),
                                  range(-vertex_neighbourhood, vertex_neighbourhood + 1, 1), indexing='ij')
            grid = numpy.c_[
                numpy.array(grid[0]).flatten(), numpy.array(grid[1]).flatten(), numpy.array(grid[2]).flatten()]
            n_grid = grid.shape[0]

            return grid, n_grid

    # Sample a volume of a specific label on a surface, by keeping only those surface vertices, the nearest voxel of
    # which is of the given label (+ of possibly additional target labels, such as white matter)
    # Allow optionally for vertices within a given voxel distance vn from the target voxels
    def sample_vol_on_surf(self, surf_path, vol_path, annot_path, out_surf_path, cras_path, ctx=None,
                           vertex_neighbourhood=1, add_lbl=[],
                           lut_path=os.path.join(os.environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt')):
        # Read the inputs
        surface = self.surface_parser.read(surf_path, False)

        annotation = self.annotation_service.annotationIO.read(annot_path)
        labels = self.annotation_service.annot_names_to_labels(annotation.region_names, ctx, lut_path)

        volume_parser = VolumeIO()
        volume = volume_parser.read(vol_path)
        ras2vox_affine_matrix = numpy.linalg.inv(volume.affine_matrix)

        cras = numpy.loadtxt(cras_path)

        grid, n_grid = self.__prepare_grid(vertex_neighbourhood)

        # Initialize the output mask:
        verts_out_mask = numpy.repeat([False], surface.vertices.shape[0])

        for label_index in range(len(labels)):
            if isinstance(ctx, basestring):
                print 'ctx-' + ctx + '-' + annotation.region_names[label_index]
            else:
                print annotation.region_names[label_index]

            # Add any additional labels
            all_labels = [labels[label_index]] + add_lbl

            # Get the indexes of the vertices corresponding to this label:
            verts_indices_of_label, = numpy.where(annotation.region_mapping[:] == label_index)
            verts_indices_of_label_size = verts_indices_of_label.size
            if verts_indices_of_label_size == 0:
                continue

            # get the vertices for current label and add cras to take them to scanner ras
            verts_of_label = surface.vertices[verts_indices_of_label, :]
            verts_of_label += numpy.repeat(numpy.expand_dims(cras, 1).T, verts_indices_of_label_size, axis=0)

            # Compute the nearest voxel coordinates using the affine transform
            ijk = numpy.round(
                ras2vox_affine_matrix.dot(numpy.c_[verts_of_label, numpy.ones(verts_indices_of_label_size)].T)[:3].T)\
                .astype('i')

            # Get the labels of these voxels:
            surf_vxls = volume.data[ijk[:, 0], ijk[:, 1], ijk[:, 2]]

            # Vertex mask to keep: those that correspond to voxels of one of the target labels
            verts_keep, = numpy.where(numpy.in1d(surf_vxls, all_labels))  # surf_vxls==lbl if only one target label
            verts_out_mask[verts_indices_of_label[verts_keep]] = True

            if vertex_neighbourhood > 0:
                # These are now the remaining indexes to be checked for neighboring voxels
                verts_indices_of_label = numpy.delete(verts_indices_of_label, verts_keep)
                ijk = numpy.delete(ijk, verts_keep, axis=0)

                for vertex_index in range(verts_indices_of_label.size):
                    # Generate the specific grid centered at the voxel ijk
                    ijk_grid = grid + numpy.tile(ijk[vertex_index, :], (n_grid, 1))

                    # Remove voxels outside the volume
                    indexes_within_limits = numpy.all([(ijk_grid[:, 0] >= 0), (ijk_grid[:, 0] < volume.dimensions[0]),
                                                       (ijk_grid[:, 1] >= 0), (ijk_grid[:, 1] < volume.dimensions[1]),
                                                       (ijk_grid[:, 2] >= 0), (ijk_grid[:, 2] < volume.dimensions[2])],
                                                      axis=0)
                    ijk_grid = ijk_grid[indexes_within_limits, :]
                    surf_vxls = volume.data[ijk_grid[:, 0], ijk_grid[:, 1], ijk_grid[:, 2]]

                    # If any of the neighbors is of the target labels include the current vertex
                    if numpy.any(numpy.in1d(surf_vxls, all_labels)):  # surf_vxls==lbl if only one target label
                        verts_out_mask[verts_indices_of_label[vertex_index]] = True

            # Vertex indexes and vertices to keep:
            verts_out_indices, = numpy.where(verts_out_mask)
            verts_out = surface.vertices[verts_out_indices]

            # TODO maybe: make sure that all voxels of this label correspond to at least one vertex.
            # Create a similar mask for faces by picking only triangles of which all 3 vertices are included
            face_out_mask = numpy.c_[
                verts_out_mask[surface.triangles[:, 0]], verts_out_mask[surface.triangles[:, 1]], verts_out_mask[
                    surface.triangles[:, 2]]].all(axis=1)
            faces_out = surface.triangles[face_out_mask]

            # The old vertices' indexes of faces have to be transformed to the new vrtx_out_inds:
            for iF in range(faces_out.shape[0]):
                for vertex_index in range(3):
                    faces_out[iF, vertex_index], = numpy.where(faces_out[iF, vertex_index] == verts_out_indices)

            # Write the output surfaces and annotations to files. Also write files with the indexes of vertices to keep.
            write_geometry(out_surf_path, verts_out, faces_out, create_stamp=None, volume_info=surface.generic_metadata)

            lab_out = annotation.region_mapping[verts_out_indices]
            self.annotation_service.annotationIO.write(out_surf_path + ".annot", lab_out, annotation.regions_color_table, annotation.region_names)

            numpy.save(out_surf_path + "-idx.npy", verts_out_indices)
            numpy.savetxt(out_surf_path + "-idx.txt", verts_out_indices, fmt='%d')
