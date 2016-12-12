# -*- coding: utf-8 -*-

import os
import numpy
import scipy
import gdist
import nibabel
# TODO maybe use classes to get an annot
from nibabel.freesurfer.io import read_geometry, write_geometry, read_annot, write_annot
from sklearn.metrics.pairwise import paired_distances
from bnm.recon.algo.service.annotation import AnnotationService


class SurfaceService(object):
    def __init__(self):
        self.annotationService = AnnotationService()

    def read_surf(hemi, name):
        surf_fname = '%s.%s' % (hemi, name)
        surf_path = os.path.join(os.environ['SUBJECTS_DIR'], os.environ['SUBJECT'], 'surf', surf_fname)
        return read_geometry(surf_path)

    def tri_area(tri):
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
            v, f = read_geometry(surf_path)
            mat_path = '%s/%s/surf/%sh.%s.gdist.mat' % (subjects_dir, subject, h, surf_name)
            mat = gdist.local_gdist_matrix(v, f.astype('<i4'), max_distance=40.0)
            scipy.io.savemat(mat_path, {'gdist': mat})

    def extract_subsurf(self, verts, faces, verts_mask):
        # These are the faces to keep...
        verts_out = verts[verts_mask, :]
        # These are the faces to keep...
        face_mask = numpy.c_[verts_mask[faces[:, 0]], verts_mask[faces[:, 1]], verts_mask[faces[:, 2]]].all(axis=1)
        faces_out = faces[face_mask]
        # ...but the old vertices' indexes of faces have to be transformed to the new verts_out_inds:
        verts_out_inds, = numpy.where(verts_mask)
        for iF in range(faces_out.shape[0]):
            for iV in range(3):
                faces_out[iF, iV], = numpy.where(faces_out[iF, iV] == verts_out_inds)
        return (verts_out, faces_out)

    def extract_mri_vol2subsurf(self, surf_path, annot_path, vol2surf_path, out_surf_path=None, out_annot_path=None,
                                ctx=None, labels=None,
                                lut_path=os.path.join(os.environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt')):
        (verts, faces, volume_info) = read_geometry(surf_path, read_metadata=True)
        vol2surf = nibabel.load(vol2surf_path)
        vol2surf = numpy.round(numpy.squeeze(vol2surf.get_data())).astype('i')
        lab, ctab, names = read_annot(annot_path)
        if labels is None:
            labels = numpy.array(self.annotationService.annot_names_to_labels(names, ctx=ctx, lut_path=lut_path))
        else:
            labels = numpy.array(labels.split()).astype('i')
        verts_mask = (v2s in labels for v2s in vol2surf)
        (verts_out, faces_out) = self.extract_subsurf(verts, faces, verts_mask)
        if os.path.exists(str(out_surf_path)):
            read_geometry(out_surf_path, verts_out, faces_out, volume_info=volume_info)
        if os.path.exists(str(out_annot_path)):
            lab = lab[verts_mask]
            write_annot(out_annot_path, lab, ctab, names)
        return (verts_out, faces_out)

    # Concatenate surfaces of specific labels to create a single annotated surface
    def aseg_surf_conc_annot(self, surf_path, out_surf_path, annot_path, labels,
                             lut_path=os.path.join(os.environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt')):
        names, ctab = self.annotationService.lut_to_annot_names_ctab(lut_path=lut_path, labels=labels)
        labels = numpy.array(labels.split()).astype('i')
        out_verts = []
        out_faces = []
        lab = []
        iL = -1
        nVerts = 0
        names_out = []
        ctab_out = []
        for lbl in labels:
            l = int(lbl)
            this_surf_path = (surf_path + "-%06d" % (l))
            if os.path.exists(this_surf_path):
                indL, = numpy.where(labels == lbl)
                names_out.append(names[indL])
                ctab_out.append(ctab[indL, :])
                iL += 1
                (verts, faces, volume_info) = read_geometry(this_surf_path, read_metadata=True)
                faces = faces + nVerts  # Update vertices indexes
                nVerts += verts.shape[0]
                out_verts.append(verts)
                out_faces.append(faces)
                lab.append(iL * numpy.ones((verts.shape[0],), dtype='int64'))
        ctab_out = numpy.squeeze(numpy.array(ctab_out).astype('i'))
        out_verts = numpy.vstack(out_verts)
        out_faces = numpy.vstack(out_faces)
        lab = numpy.hstack(lab)
        write_geometry(out_surf_path, out_verts, out_faces, create_stamp=None, volume_info=volume_info)
        write_annot(annot_path, lab, ctab_out, names_out)

    # It returns a sparse matrix of the connectivity among the vertices of a surface
    # mode: "sparse" (default) or "2D"
    def vertex_connectivity(self, v, f, mode="sparse", metric=None):
        # Get all pairs of vertex indexes that appear in each face
        f = numpy.r_[f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]
        # Remove repetitions
        f = numpy.vstack(set(map(tuple, f)))
        # Mark all existing pairs to 1
        nV = v.shape[0]
        nF = f.shape[0]
        from scipy.sparse import csr_matrix
        if metric is None:
            con = csr_matrix((numpy.ones((nF,)), (f[:, 0], f[:, 1])), shape=(nV, nV))
            if mode != "sparse":
                # Create non-sparse matrix
                con = con.todense()
        else:
            d = paired_distances(v[f[:, 0]], v[f[:, 1]], metric)
            if mode == "sparse":
                # Create sparse matrix
                con = csr_matrix((d, (f[:, 0], f[:, 1])), shape=(nV, nV))
        return con

    # Sample a volume of a specific label on a surface, by keeping
    # only those surface vertices, the nearest voxel of which is of the given label (+ of possibly additional target labels, such as white matter)
    # Allow optionally for vertices within a given voxel distance vn from the target voxels
    def sample_vol_on_surf(self, surf_path, vol_path, annot_path, out_surf_path, cras_path, ctx=None, vn=1, add_lbl=[],
                           lut_path=os.path.join(os.environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt')):
        # Read the surface...
        (verts, faces, volume_info) = read_geometry(surf_path, read_metadata=True)
        # ...and its annotation:
        lab, ctab, names = read_annot(annot_path)
        # Get the region names of these labels:
        labels = self.annotationService.annot_names_to_labels(names, ctx, lut_path)
        nLbl = len(labels)
        # Read the volume...
        volume = nibabel.load(vol_path)
        # ...and get its data
        vol = volume.get_data()
        vol_shape = vol.shape
        # ...and invert its vox2ras transform
        vox2ras = volume.affine
        xyz2ijk = numpy.linalg.inv(vox2ras)
        # Read the cras
        cras = numpy.loadtxt(cras_path)
        # Prepare grid if needed for possible use:
        if vn > 0:
            grid = numpy.meshgrid(range(-vn, vn + 1, 1), range(-vn, vn + 1, 1), range(-vn, vn + 1, 1), indexing='ij')
            grid = numpy.c_[
                numpy.array(grid[0]).flatten(), numpy.array(grid[1]).flatten(), numpy.array(grid[2]).flatten()]
            nGrid = grid.shape[0]
        # Initialize the output mask:
        verts_out_mask = numpy.repeat(False, verts.shape[0])
        for iL in range(nLbl):
            if isinstance(ctx, basestring):
                print 'ctx-' + ctx + '-' + names[iL]
            else:
                print names[iL]
            # Form the target labels by adding to the input label list any additional labels, if any
            lbl = [labels[iL]] + add_lbl
            # Get the indexes of the vertices of this label:
            verts_lbl_inds, = numpy.where(lab[:] == iL)
            nVlbl = verts_lbl_inds.size
            if nVlbl == 0:
                continue
            # Apply the affine transform to the selected surface vertices, and get integer indices
            # of the corresponding nearest voxels in the vol
            # Get the specific vertices in tkras coordinates...
            verts_lbl = verts[verts_lbl_inds, :]
            # ...add the cras to take them to scanner ras...
            verts_lbl += numpy.repeat(numpy.expand_dims(cras, 1).T, nVlbl, axis=0)
            # ...and compute the nearest voxel coordinates
            ijk = numpy.round(xyz2ijk.dot(numpy.c_[verts_lbl, numpy.ones(nVlbl)].T)[:3].T).astype('i')
            # Get the labels of these voxels:
            surf_vxls = vol[ijk[:, 0], ijk[:, 1], ijk[:, 2]]
            # Vertex mask to keep: those that correspond to voxels of one of the target labels
            verts_keep, = numpy.where(numpy.in1d(surf_vxls, lbl))  # surf_vxls==lbl if only one target label
            verts_out_mask[verts_lbl_inds[verts_keep]] = True
            if vn > 0:
                # These are now the remaining indexes to be checked for neighboring voxels
                verts_lbl_inds = numpy.delete(verts_lbl_inds, verts_keep)
                ijk = numpy.delete(ijk, verts_keep, axis=0)
                for iV in range(verts_lbl_inds.size):
                    # Generate the specific grid centered at the voxel ijk
                    ijk_grid = grid + numpy.tile(ijk[iV, :], (nGrid, 1))
                    # Remove voxels outside the volume
                    indexes_within_limits = numpy.all([(ijk_grid[:, 0] >= 0), (ijk_grid[:, 0] < vol_shape[0]),
                                                       (ijk_grid[:, 1] >= 0), (ijk_grid[:, 1] < vol_shape[1]),
                                                       (ijk_grid[:, 2] >= 0), (ijk_grid[:, 2] < vol_shape[2])],
                                                      axis=0)
                    ijk_grid = ijk_grid[indexes_within_limits, :]
                    # Get the labels of these voxels:
                    surf_vxls = vol[ijk_grid[:, 0], ijk_grid[:, 1], ijk_grid[:, 2]]
                    # If any of the neighbors is of the target labels...
                    if numpy.any(numpy.in1d(surf_vxls, lbl)):  # surf_vxls==lbl if only one target label
                        # ...include this vertex
                        verts_out_mask[verts_lbl_inds[iV]] = True
            # Vertex indexes to keep:
            verts_out_inds, = numpy.where(verts_out_mask)
            # These are the vertices to keep
            verts_out = verts[verts_out_inds]
            # TODO maybe: make sure that all voxels of this label correspond to at least one vertex.
            # Create a similar mask for faces by picking only triangles
            # of which all 3 vertices are included
            face_out_mask = numpy.c_[
                verts_out_mask[faces[:, 0]], verts_out_mask[faces[:, 1]], verts_out_mask[faces[:, 2]]].all(axis=1)
            # These are the faces to keep...
            faces_out = faces[face_out_mask]
            # ...but the old vertices' indexes of faces have to be transformed to the new vrtx_out_inds:
            for iF in range(faces_out.shape[0]):
                for iV in range(3):
                    faces_out[iF, iV], = numpy.where(faces_out[iF, iV] == verts_out_inds)
            # Write the output surfaces to a file
            write_geometry(out_surf_path, verts_out, faces_out, create_stamp=None, volume_info=volume_info)
            # Create and write output annotations to files
            lab_out = lab[verts_out_inds]
            write_annot(out_surf_path + ".annot", lab_out, ctab, names)
            # Write files with the indexes of vertices to keep
            numpy.save(out_surf_path + "-idx.npy", verts_out_inds)
            numpy.savetxt(out_surf_path + "-idx.txt", verts_out_inds, fmt='%d')
