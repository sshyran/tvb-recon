# -*- coding: utf-8 -*-

import os
import nibabel
import numpy
from bnm.recon.algo.service.annotation import AnnotationService

class VolumeService(object):
    def __init__(self):
        self.annotationservice = AnnotationService()

    # Separate the voxels of the outer surface of a structure, from the inner ones
    # Default behavior: surface voxels retain their label, inner voxels get the label 0,
    # and the input file is overwritten by the output
    def vol_to_ext_surf_vol(self, in_vol_path, labels=None, hemi=None, out_vol_path=None, labels_surf=None,
                            labels_inner='0'):
        # Set the target labels:
        (labels, nLbl) = self.annotationservice.read_input_labels(labels=labels, hemi=hemi)
        # Set the labels for the surfaces
        if labels_surf is None:
            labels_surf = labels
        else:
            # Read the surface labels
            labels_surf = numpy.array(labels_surf.split()).astype('i')
            # ...and make sure there is one for each label
            if len(labels_surf) == 1:
                labels_surf = numpy.repeat(labels_inner, nLbl).tolist()
            elif len(labels_surf) != nLbl:
                print "Output labels for surface voxels are neither of length 1 nor of length equal to the one of target labels"
            else:
                labels_surf = labels_surf.tolist()
        # Read the inner, non-surface labels
        labels_inner = numpy.array(labels_inner.split()).astype('i')
        # ...and make sure there is one for each label
        if len(labels_inner) == 1:
            labels_inner = numpy.repeat(labels_inner, nLbl).tolist()
        elif len(labels_inner) != nLbl:
            print "Output labels for inner voxels are neither of length 1 nor of length equal to the one of the target labels"
        else:
            labels_inner = labels_inner.tolist()
        # Read the input volume...
        volume = nibabel.load(in_vol_path)
        # ...and get its data
        vol = volume.get_data()
        vol_shape = vol.shape
        # Neigbors' grid sharing a face
        borderGrid = numpy.c_[numpy.identity(3), -numpy.identity(3)].T.astype('i')
        nBorder = 6
        # Initialize output volume array
        out_vol = numpy.array(vol)
        # Initialize output indexes
        out_ijk = []
        # For each target label:
        for iL in range(nLbl):
            # this label
            lbl = labels[iL]
            # Get the indexes of all voxels of this label:
            ii, jj, kk = numpy.where(vol == lbl)
            # and for each voxel
            for iV in range(ii.size):
                # indexes of this voxel:
                (i, j, k) = (ii[iV], jj[iV], kk[iV])
                # Create the neighbors' grid sharing a face
                ijk_grid = borderGrid + numpy.tile(numpy.array([i, j, k]), (nBorder, 1))
                # Remove voxels outside the image
                inds_inside_image = numpy.all([(ijk_grid[:, 0] >= 0), (ijk_grid[:, 0] < vol_shape[0]),
                                               (ijk_grid[:, 1] >= 0), (ijk_grid[:, 1] < vol_shape[1]),
                                               (ijk_grid[:, 2] >= 0), (ijk_grid[:, 2] < vol_shape[2])],
                                              axis=0)
                ijk_grid = ijk_grid[inds_inside_image, :]
                try:
                    # If all face neighbors are of the same label...
                    if numpy.all(vol[ijk_grid[:, 0], ijk_grid[:, 1], ijk_grid[:, 2]] == numpy.tile(vol[i, j, k],
                                                                                                   (nBorder, 1))):
                        # ...set this voxel to the corresponding inner target label
                        out_vol[i, j, k] = labels_inner[iL]
                    else:
                        # ...set this voxel to the corresponding surface target label
                        out_vol[i, j, k] = labels_surf[iL]
                        out_ijk.append([i, j, k])
                except ValueError:  # empty grid
                    print "Error at voxel (" + str(i) + "," + str(j) + "," + str(k) + ") of label " + str(lbl) + ":"
                    print "It appears to have no common-face neighbors inside the image!"
                    return
                    # Create the new volume and save it
        out_volume = nibabel.Nifti1Image(out_vol, volume.affine, header=volume.header)
        if out_vol_path == None:
            # Overwrite volume
            out_vol_path = in_vol_path
            # Save a new volume
        nibabel.save(out_volume, out_vol_path)
        # ...and the output indexes that survived masking
        out_ijk = numpy.vstack(out_ijk)
        filepath = os.path.splitext(out_vol_path)[0]
        numpy.save(filepath + "-idx.npy", out_ijk)
        numpy.savetxt(filepath + "-idx.txt", out_ijk, fmt='%d')

    # Identify the voxels that our neighbors with a voxel distance vn, to a mask volume,
    # with a mask threshold of th
    # Default behavior: we assume a binarized mask and set th=0.999,
    # no neigbhors search, only looking at the exact voxel position, i.e., vn=0.
    # and mask voxels retain their label, no mask voxels get a label of 0
    def mask_to_vol(self, in_vol_path, mask_vol_path, out_vol_path=None, labels=None, hemi=None, vol2mask_path=None, vn=1,
                    th=0.999, labels_mask=None, labels_nomask='0'):
        # Set the target labels:
        (labels, nLbl) = self.annotationservice.read_input_labels(labels=labels, hemi=hemi)
        # Set the labels for the selected voxels
        if labels_mask is None:
            labels_mask = labels
        else:
            # Read the labels
            labels_mask = numpy.array(labels_mask.split()).astype('i')
            # ...and make sure there is one for each label
            if len(labels_mask) == 1:
                labels_mask = numpy.repeat(labels_mask, nLbl).tolist()
            elif len(labels_mask) != nLbl:
                print "Output labels for selected voxels are neither of length 1 nor of length equal to the one of target labels"
            else:
                labels_mask = labels_mask.tolist()
        # Read the excluded labels
        labels_nomask = numpy.array(labels_nomask.split()).astype('i')
        # ...and make sure there is one for each label
        if len(labels_nomask) == 1:
            labels_nomask = numpy.repeat(labels_nomask, nLbl).tolist()
        elif len(labels_nomask) != nLbl:
            print "Output labels for excluded voxels are neither of length 1 nor of length equal to the one of the target labels"
        else:
            labels_nomask = labels_nomask.tolist()
        # Read the target volume...
        volume = nibabel.load(in_vol_path)
        # ...and get its data
        vol = volume.get_data()
        # ...and its affine transform
        # ijk2xyz_vol = volume.affine
        # Read the mask volume...
        mask_vol = nibabel.load(mask_vol_path)
        # ...and get its data
        mask = mask_vol.get_data()
        mask_shape = mask.shape
        # ...and invert its affine transform
        # xyz2ijk_mask = numpy.linalg.inv(mask_vol.affine)
        # Finally compute the transform from vol ijk to mask ijk:
        ijk2ijk = numpy.identity(4)
        # If vol and mask are not in the same space:
        if os.path.exists(str(vol2mask_path)):
            # read the xyz2xyz transform...
            xyz2xyz = numpy.loadtxt(vol2mask_path)
            # ...and apply it to the inverse mask affine transform to get an ijk2ijk transform:
            ijk2ijk = volume.affine.dot(numpy.dot(xyz2xyz, numpy.linalg.inv(mask_vol.affine)))
            # Construct a grid template of voxels +/- vn voxels around each ijk voxel,
        # sharing at least a corner
        grid = numpy.meshgrid(range(-vn, vn + 1, 1), range(-vn, vn + 1, 1), range(-vn, vn + 1, 1), indexing='ij')
        grid = numpy.c_[numpy.array(grid[0]).flatten(), numpy.array(grid[1]).flatten(), numpy.array(grid[2]).flatten()]
        nGrid = grid.shape[0]
        # Initialize the output volume
        out_vol = numpy.array(vol)
        # Initialize output indexes
        out_ijk = []
        # For each target label:
        for iL in range(nLbl):
            lbl = labels[iL]
            # Get the indexes of all voxels of this label:
            ii, jj, kk = numpy.where(vol == lbl)
            # and for each voxel
            for iV in range(ii.size):
                # indexes of this voxel:
                (i, j, k) = (ii[iV], jj[iV], kk[iV])
                # TODO if necessary: deal with voxels at the edge of the image, such as brain stem ones...
                #           #Check if it is a border voxel:
                #           if any([(i==0), (i==mask_shape[0]-1),(j==0), (j==mask_shape[0]-1),(k==0), (k==mask_shape[0]-1)]):
                #               #set this voxel to the 0 label
                #               mask_shape[i,j,k]=0
                #               continue
                # ...get the corresponding voxel in the mask volume:
                ijk = numpy.round(ijk2ijk.dot(numpy.array([i, j, k, 1]))[:3]).astype('i')
                # Make sure this point is within image limits
                for cc in range(3):
                    if ijk[cc] < 0:
                        ijk[cc] = 0
                    elif ijk[cc] >= mask_shape[cc]:
                        ijk[cc] = mask_shape[cc] - 1
                # If this is a voxel to keep, set it so...
                if (mask[ijk[0], ijk[1], ijk[2]] >= th):
                    out_vol[i, j, k] = labels_mask[iL]
                    out_ijk.append([i, j, k])
                elif vn > 0:
                    # ...if not, and as long as vn>0...
                    # ...check whether any of its vn neighbors is a mask voxel
                    # Generate the specific grid centered at the vertex ijk
                    ijk_grid = grid + numpy.tile(ijk, (nGrid, 1))
                    # Remove voxels outside the mask volume
                    indexes_within_limits = numpy.all([(ijk_grid[:, 0] >= 0), (ijk_grid[:, 0] < mask_shape[0]),
                                                       (ijk_grid[:, 1] >= 0), (ijk_grid[:, 1] < mask_shape[1]),
                                                       (ijk_grid[:, 2] >= 0), (ijk_grid[:, 2] < mask_shape[2])],
                                                      axis=0)
                    ijk_grid = ijk_grid[indexes_within_limits, :]
                    try:
                        # If none of these points is a mask point:
                        if (mask[ijk_grid[:, 0], ijk_grid[:, 1], ijk_grid[:, 2]] < th).all():
                            out_vol[i, j, k] = labels_nomask[iL]
                        else:  # if any of them is a mask point:
                            out_vol[i, j, k] = labels_mask[iL]
                            out_ijk.append([i, j, k])
                    except ValueError:  # empty grid
                        print "Error at voxel (" + str(i) + "," + str(j) + "," + str(k) + "):"
                        print "It appears to have no common-face neighbors inside the image!"
                        return
                else:
                    out_vol[i, j, k] = labels_nomask[iL]
        # Create the new volume and save it
        out_volume = nibabel.Nifti1Image(out_vol, volume.affine, header=volume.header)
        if out_vol_path == None:
            # Overwrite volume
            out_vol_path = in_vol_path
        nibabel.save(out_volume, out_vol_path)
        # ...and the output indexes that survived masking
        out_ijk = numpy.vstack(out_ijk)
        filepath = os.path.splitext(out_vol_path)[0]
        numpy.save(filepath + "-idx.npy", out_ijk)
        numpy.savetxt(filepath + "-idx.txt", out_ijk, fmt='%d')

    def label_with_dilation(self, to_label_nii_fname, dilated_nii_fname, out_nii_fname):
        "Label one nifti with its dilation, cf seeg-ct.sh"
        # TODO could make dilation with ndimage also.
        import scipy.ndimage
        mask = nibabel.load(to_label_nii_fname)
        dil_mask = nibabel.load(dilated_nii_fname)
        lab, n = scipy.ndimage.label(dil_mask.get_data())
        print('[label_with_dilation] %d objects found.' % (n,))
        lab_mask_nii = nibabel.nifti1.Nifti1Image(lab * mask.get_data(), mask.affine)
        nibabel.save(lab_mask_nii, out_nii_fname)

    def label_vol_from_tdi(self, tdi_nii_fname, out_fname, lo=0.5):
        "Make label volume from tckmap output."
        # Load tdi_ends volume:
        nii = nibabel.load(tdi_nii_fname)
        # Copy its data...
        tdi = nii.get_data().copy()
        # and mask them to get the voxels of tract ends
        mask = tdi > lo
        # (all other voxels ->0)
        tdi[~mask] = 0
        # Assign them with integer labels starting from 1
        tdi[mask] = numpy.r_[1:mask.sum() + 1]
        # Write tdi_lbl to file
        out_nii = nibabel.nifti1.Nifti1Image(tdi, nii.affine)
        nibabel.save(out_nii, out_fname)

    # It removes network nodes with zero connectivity, and returns a symmetric connectivity matrix
    # Inputs:
    #    - the tdi_lbl.nii volume path
    #    -a .csv file path, output of Mrtrix3 tck2connectome
    # Outputs:
    #   - the symmetric no-zero-connections connectivity matrix saved as .npy
    #   - the tdi_lbl.nii volume with the removed voxel nodes to 0 and the labels
    #    updated
    # Optionally: if the tract length matrix is in the input, it is also processed
    def remove_zero_connectivity_nodes(self, node_vol_path, con_mat_path, tract_length_path=None):
        # Read input files:
        # Nodes' volume (nii):
        node_vol = nibabel.load(node_vol_path)
        vol = node_vol.get_data()
        # Connectivity matrix (.csv)
        con = numpy.array(numpy.genfromtxt(con_mat_path, dtype='int64'))
        # Make it symmetric:
        con = con + con.T
        # Sum rows to get the total connectivity per node
        consum = numpy.sum(con, axis=0)
        # Index of nodes to keep
        ii = consum > 0
        # Select only the specified columns and rows from the connectivity matrix:
        con = con[ii, :][:, ii]
        # Write output files
        numpy.save(os.path.splitext(con_mat_path)[0] + ".npy", con)
        numpy.savetxt(con_mat_path, con)
        # If there is also a tract length file:
        if os.path.exists(str(tract_length_path)):
            # read it:
            con = numpy.array(numpy.genfromtxt(tract_length_path, dtype='int64'))
            # Select only the specified columns and rows from the connectivity matrix:
            con = con[ii, :][:, ii]
            # Write output files
            numpy.save(os.path.splitext(tract_length_path)[0] + ".npy", con)
            numpy.savetxt(tract_length_path, con)
        else:
            print tract_length_path + " is not a valid path"
        # Index of nodes to remove
        ii, = numpy.where(~ii)
        ii = ii + 1
        nKeep = con.shape[0]
        # Remove:
        for iR in ii:
            vol[vol == iR] = 0
        # Update remaining indexes
        vol[vol > 0] = numpy.r_[1:(nKeep + 1)]
        # Write the updated volume file:
        node_vol = nibabel.Nifti1Image(vol, node_vol.affine, header=node_vol.header)
        nibabel.save(node_vol, node_vol_path)

    # It receives a binary connectivity matrix, and outputs a node connectivity
    # similarity or distance  matrix
    # con_mat_path: path to connectivity file
    # metric: default "cosine"
    # mode: "sim" or "dist" for similarity or distance output
    def node_connectivity_metric(self, con_mat_path, metric="cosine", mode='sim', out_consim_path=None):
        from scipy.spatial.distance import pdist, squareform
        # Read iput file
        con = numpy.load(con_mat_path)
        # Calculate distance metric
        con = squareform(pdist(con, metric=metric))
        # If similarity is required,...
        if mode == 'sim':
            # ...calculate it:
            con = 1 - con
        if out_consim_path is not None:
            numpy.save(out_consim_path, con)
        return con

    def simple_label_config(self, aparc_fname, out_fname):
        "Rewrite label volume to have contiguous values like mrtrix' labelconfig."
        aparc = nibabel.load(aparc_fname)
        vol = aparc.get_data()
        uval = numpy.unique(vol)
        uval_map = numpy.r_[:uval.max() + 1]
        uval_map[uval] = numpy.r_[:uval.size]
        uvol = uval_map[vol]
        uparc = nibabel.nifti1.Nifti1Image(uvol, aparc.affine)
        nibabel.save(uparc, out_fname)
