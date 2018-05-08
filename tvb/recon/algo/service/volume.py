# -*- coding: utf-8 -*-

import os
from typing import Optional, Union
import numpy
import scipy.ndimage
from tvb.recon.dax import AtlasSuffix
import nibabel
from tvb.recon.logger import get_logger
from tvb.recon.algo.service.utils import execute_command
from tvb.recon.algo.service.annotation import AnnotationService
from tvb.recon.io.factory import IOUtils
from tvb.recon.model.volume import Volume
from tvb.recon.model.constants import NPY_EXTENSION


class VolumeService(object):
    logger = get_logger(__name__)

    def __init__(self):
        self.annotation_service = AnnotationService()

    def gen_label_volume_from_labels_inds(self, values: Union[numpy.ndarray, list],
                                          input_label_volume_file: os.PathLike, output_label_volume_file: os.PathLike) \
            -> nibabel.nifti1.Nifti1Image:
        """
        This function creates a new volume matrix from an input label_volume matrix
        by setting values[i] to all voxels where label_volume == i
        """

        label_nii = nibabel.load(input_label_volume_file)
        label_volume = label_nii.get_data()

        new_volume = numpy.zeros(label_volume.shape)
        new_volume[:, :] = numpy.nan

        for i, value in enumerate(values):
            region = i + 1
            mask = label_volume[:, :, :] == region
            new_volume[mask] = value

        # TODO: I don't know what this is... I have to ask Viktor...
        def add_min_max(volume):
            """Ugly hack to deal with the MRtrix bug (?) that causes MRview to crop min/max values"""
            volume[0, 0, 0] = numpy.nanmin(volume) - 1
            volume[0, 0, 1] = numpy.nanmax(volume) + 1

        # add_min_max(new_volume)
        new_nii = nibabel.Nifti1Image(new_volume, label_nii.affine)
        nibabel.save(new_nii, output_label_volume_file)

        return new_nii

    def gen_label_volume_from_coords(self, coords: Union[os.PathLike, numpy.ndarray],
                                     ref_volume_file: os.PathLike, out_volume_file: os.PathLike,
                                     values: Union[numpy.ndarray, list]=None,
                                     labels: Union[os.PathLike, numpy.ndarray, list]=None,
                                     skip_missing: bool=False, dist: int=0) \
            -> numpy.ndarray:
        """
        # Create and save a new nifti label volume of similar shape to a reference volume
        # by setting input values at input positions (optionally + & - (int) dist)
        # after applying to them the ref_aff transform

        :param values: vector of values to be included in the output volume, numpy.ndarray of shape (n_vals, )
        :param coords: array or list of 3D point coordinates, numpy.ndarray of shape (n_coords, 3), or file to read array
        :param labels: array or list of labels'names, or file to read them
        :param ref_volume_file: path to nifti volume
        :param out_volume_file: file path for the output nifti volume to be written
        :param skip_missing: flag
        :param dist: integer indicating the size of the neighborhood around the coords' positions to be labeled
        :return: output nifti volume
        """
        ref_volume = nibabel.load(ref_volume_file)

        if os.path.isfile(str(coords)):
            coords = numpy.genfromtxt(coords, dtype=float, usecols=(1, 2, 3))
        n_coords = coords.shape[0]

        if labels is None:
            labels = range(n_coords)
        else:
            if os.path.isfile(str(labels)):
                labels = list(numpy.genfromtxt(labels, dtype=str, usecols=(0,)))
            else:
                labels = list(labels)

        if values is None:
            values = numpy.ones((n_coords, ))

        positions = numpy.zeros((len(values), 3))
        for i, label in enumerate(labels):
            try:
                contact_ind = labels.index(label)
                coord = coords[contact_ind, :]
            except ValueError:
                coord = numpy.array([numpy.nan, numpy.nan, numpy.nan])

            positions[i, :] = coord

        missing_mask = numpy.isnan(positions[:, 0])
        if skip_missing:
            values = values[~missing_mask]
            positions = positions[~missing_mask, :]
        else:
            if numpy.any(missing_mask):
                raise ValueError("Missing contact position(s) for: %s." % ", ".join([
                    name for name, missing in zip(labels, missing_mask) if missing]))

            # numpy.array(names)[missing_mask]))

        new_volume = numpy.zeros(ref_volume.shape)
        # TODO: Find out the use of the following commented line:
        # new_volume[:, :] = numpy.nan

        kx, ky, kz = numpy.mgrid[-dist:dist + 1, -dist:dist + 1, -dist:dist + 1]

        for val, pos in zip(values, positions):
            ix, iy, iz = numpy.linalg.solve(ref_volume.affine, numpy.append(pos, 1.0))[0:3].astype(int)
            # new_volume[inds[0], inds[1], inds[2]] = val
            new_volume[ix + kx, iy + ky, iz + kz] = val

        # add_min_max(new_volume)

        new_nii = nibabel.Nifti1Image(new_volume, ref_volume.affine)
        nibabel.save(new_nii, out_volume_file)

        return nibabel.nifti1.Nifti1Image

    def vol_to_ext_surf_vol(self, in_vol_path: os.PathLike, labels: Optional[Union[numpy.ndarray, list]]=None,
                            ctx: Optional[os.PathLike]=None, out_vol_path: Optional[os.PathLike]=None,
                            labels_surf: Optional[Union[numpy.ndarray, list]]=None, labels_inner: str='0'):
        """
        Separate the voxels of the outer surface of a structure, from the inner ones. Default behavior: surface voxels
        retain their label, inner voxels get the label 0, and the input file is overwritten by the output.
        """

        labels = self.annotation_service.read_input_labels(
            labels=labels, ctx=ctx)
        number_of_labels = len(labels)
        # Set the labels for the surfaces
        if labels_surf is None:
            labels_surf = labels
        else:
            # Read the surface labels and make sure there is one for each label
            labels_surf = numpy.array(labels_surf.split()).astype('i')
            if len(labels_surf) == 1:
                labels_surf = numpy.repeat(
                    labels_inner, number_of_labels).tolist()
            elif len(labels_surf) != number_of_labels:
                self.logger.warning(
                    "Output labels for surface voxels are neither of length "
                    "1 nor of length equal to the one of target labels.")
                return
            else:
                labels_surf = labels_surf.tolist()
        # Read the inner, non-surface labels
        labels_inner = numpy.array(labels_inner.split()).astype('i')
        # ...and make sure there is one for each label
        if len(labels_inner) == 1:
            labels_inner = numpy.repeat(
                labels_inner, number_of_labels).tolist()
        elif len(labels_inner) != number_of_labels:
            self.logger.warning(
                "Output labels for inner voxels are neither of length 1 nor "
                "of length equal to the one of the target labels.")
            return
        else:
            labels_inner = labels_inner.tolist()

        # Read the input volume...
        volume = IOUtils.read_volume(in_vol_path)

        # Neigbors' grid sharing a face
        eye3 = numpy.identity(3)
        border_grid = numpy.c_[eye3, -eye3].T.astype('i')
        n_border = 6

        out_volume = Volume(numpy.array(volume.data),
                            volume.affine_matrix, volume.header)

        # Initialize output indexes
        out_ijk = []

        for label_index in range(number_of_labels):
            current_label = labels[label_index]
            # Get the indexes of all voxels of this label:
            label_volxels_i, label_voxels_j, label_voxels_k = numpy.where(
                volume.data == current_label)
            # and for each voxel
            for voxel_index in range(label_volxels_i.size):
                # indexes of this voxel:
                current_voxel_i, current_voxel_j, current_voxel_k = \
                    label_volxels_i[voxel_index], label_voxels_j[
                        voxel_index], label_voxels_k[voxel_index]
                # Create the neighbors' grid sharing a face
                ijk_grid = border_grid + \
                           numpy.tile(numpy.array(
                               [current_voxel_i, current_voxel_j, current_voxel_k]), (n_border, 1))
                # Remove voxels outside the image
                indices_inside_image = numpy.all([(ijk_grid[:, 0] >= 0), (ijk_grid[:, 0] < volume.dimensions[0]),
                                                  (ijk_grid[:, 1] >= 0), (ijk_grid[
                                                                          :, 1] < volume.dimensions[1]),
                                                  (ijk_grid[:, 2] >= 0), (ijk_grid[:, 2] < volume.dimensions[2])],
                                                 axis=0)
                ijk_grid = ijk_grid[indices_inside_image, :]
                try:
                    # If all face neighbors are of the same label...
                    if numpy.all(volume.data[ijk_grid[:, 0], ijk_grid[:, 1], ijk_grid[:, 2]] == numpy.tile(
                            volume.data[current_voxel_i,
                                        current_voxel_j, current_voxel_k],
                            (n_border, 1))):
                        # ...set this voxel to the corresponding inner target label
                        out_volume.data[current_voxel_i, current_voxel_j,
                                        current_voxel_k] = labels_inner[label_index]
                    else:
                        # ...set this voxel to the corresponding surface target label
                        out_volume.data[current_voxel_i, current_voxel_j,
                                        current_voxel_k] = labels_surf[label_index]
                        out_ijk.append(
                            [current_voxel_i, current_voxel_j, current_voxel_k])
                except ValueError:  # empty grid
                    self.logger.error("Error at voxel ( %s, %s, %s ) of label %s: It appears to have no common-face "
                                      "neighbors inside the image!", str(
                        current_voxel_i), str(current_voxel_j),
                                      str(current_voxel_k), str(current_label))
                    return

        if out_vol_path is None:
            out_vol_path = in_vol_path

        IOUtils.write_volume(out_vol_path, out_volume)

        # save the output indexes that survived masking
        out_ijk = numpy.vstack(out_ijk)
        filepath = os.path.splitext(out_vol_path)[0]
        numpy.save(filepath + "-idx.npy", out_ijk)
        numpy.savetxt(filepath + "-idx.txt", out_ijk, fmt='%d')

    def mask_to_vol(self, in_vol_path: os.PathLike, mask_vol_path: os.PathLike,
                    out_vol_path: Optional[os.PathLike]=None, labels: Optional[Union[numpy.ndarray, list]]=None,
                    ctx: Optional[str]=None, vol2mask_path: Optional[os.PathLike]=None, vn: int=1, th: float=0.999,
                    labels_mask: Optional[os.PathLike]=None, labels_nomask: str='0'):
        """
        Identify the voxels that are neighbors within a voxel distance vn, to a mask volume, with a mask threshold of th
        Default behavior: we assume a binarized mask and set th=0.999, no neighbors search, only looking at the exact
        voxel position, i.e., vn=0. Accepted voxels retain their label, whereas rejected ones get a label of 0
        """

        # Set the target labels:
        labels = self.annotation_service.read_input_labels(
            labels=labels, ctx=ctx)
        number_of_labels = len(labels)
        # Set the labels for the selected voxels
        if labels_mask is None:
            labels_mask = labels

        else:
            # Read the labels and make sure there is one for each label
            labels_mask = numpy.array(labels_mask.split()).astype('i')

            if len(labels_mask) == 1:
                labels_mask = numpy.repeat(
                    labels_mask, number_of_labels).tolist()

            elif len(labels_mask) != number_of_labels:
                self.logger.warning("Output labels for selected voxels are neither of length 1 nor of length equal to "
                                    "the one of target labels")
                return

            else:
                labels_mask = labels_mask.tolist()

        # Read the excluded labels and make sure there is one for each label
        labels_nomask = numpy.array(labels_nomask.split()).astype('i')
        if len(labels_nomask) == 1:
            labels_nomask = numpy.repeat(
                labels_nomask, number_of_labels).tolist()

        elif len(labels_nomask) != number_of_labels:
            self.logger.warning("Output labels for excluded voxels are neither of length 1 nor of length equal to the "
                                "one of the target labels")
            return

        else:
            labels_nomask = labels_nomask.tolist()

        volume = IOUtils.read_volume(in_vol_path)

        mask_vol = IOUtils.read_volume(mask_vol_path)

        # Compute the transform from vol ijk to mask ijk:
        ijk2ijk = numpy.identity(4)

        # If vol and mask are not in the same space:
        if os.path.exists(str(vol2mask_path)):
            # read the xyz2xyz transform and apply it to the inverse mask
            # affine transform to get an ijk2ijk transform.
            xyz2xyz = numpy.loadtxt(vol2mask_path)
            ijk2ijk = volume.affine_matrix.dot(
                numpy.dot(xyz2xyz, numpy.linalg.inv(mask_vol.affine_matrix)))

        # Construct a grid template of voxels +/- vn voxels around each ijk
        # voxel, sharing at least a corner
        grid = numpy.meshgrid(list(range(-vn, vn + 1, 1)), list(
            range(-vn, vn + 1, 1)), list(range(-vn, vn + 1, 1)), indexing='ij')
        grid = numpy.c_[numpy.array(grid[0]).flatten(), numpy.array(
            grid[1]).flatten(), numpy.array(grid[2]).flatten()]
        n_grid = grid.shape[0]

        out_volume = Volume(numpy.array(volume.data),
                            volume.affine_matrix, volume.header)

        # Initialize output indexes
        out_ijk = []

        # For each target label:
        for label_index in range(number_of_labels):
            current_label = labels[label_index]
            # Get the indexes of all voxels of this label:
            label_voxels_i, label_voxels_j, label_voxels_k = numpy.where(
                volume.data == current_label)

            for voxel_index in range(label_voxels_i.size):
                current_voxel_i, current_voxel_j, current_voxel_k = \
                    label_voxels_i[voxel_index], label_voxels_j[
                        voxel_index], label_voxels_k[voxel_index]
                # TODO if necessary: deal with voxels at the edge of the image, such as brain stem ones...
                #     if any([(i==0), (i==mask_shape[0]-1),(j==0), (j==mask_shape[0]-1),(k==0), (k==mask_shape[0]-1)]):
                #               mask_shape[i,j,k]=0
                #               continue

                # ...get the corresponding voxel in the mask volume:
                ijk = numpy.round(ijk2ijk.dot(numpy.array(
                    [current_voxel_i, current_voxel_j, current_voxel_k, 1]))[:3]).astype('i')

                # Make sure this point is within image limits
                for cc in range(3):
                    if ijk[cc] < 0:
                        ijk[cc] = 0

                    elif ijk[cc] >= mask_vol.dimensions[cc]:
                        ijk[cc] = mask_vol.dimensions[cc] - 1

                # If this is a voxel to keep, set it so...
                if mask_vol.data[ijk[0], ijk[1], ijk[2]] >= th:
                    out_volume.data[current_voxel_i, current_voxel_j,
                                    current_voxel_k] = labels_mask[label_index]
                    out_ijk.append(
                        [current_voxel_i, current_voxel_j, current_voxel_k])

                elif vn > 0:
                    # If not, and as long as vn>0 check whether any of its vn neighbors is a mask voxel.
                    # Generate the specific grid centered at the vertex ijk
                    ijk_grid = grid + numpy.tile(ijk, (n_grid, 1))

                    # Remove voxels outside the mask volume
                    indexes_within_limits = numpy.all([(ijk_grid[:, 0] >= 0), (ijk_grid[:, 0] < mask_vol.dimensions[0]),
                                                       (ijk_grid[:, 1] >= 0), (ijk_grid[
                                                                               :, 1] < mask_vol.dimensions[1]),
                                                       (ijk_grid[:, 2] >= 0),
                                                       (ijk_grid[:, 2] < mask_vol.dimensions[2])],
                                                      axis=0)
                    ijk_grid = ijk_grid[indexes_within_limits, :]

                    try:
                        # If none of these points is a mask point:
                        if (mask_vol.data[ijk_grid[:, 0], ijk_grid[
                                                          :, 1], ijk_grid[:, 2]] < th).all():
                            out_volume.data[
                                current_voxel_i, current_voxel_j, current_voxel_k] = labels_nomask[label_index]

                        else:  # if any of them is a mask point:
                            out_volume.data[
                                current_voxel_i, current_voxel_j, current_voxel_k] = labels_mask[label_index]
                            out_ijk.append(
                                [current_voxel_i, current_voxel_j, current_voxel_k])

                    except ValueError:  # empty grid
                        self.logger.error("Error at voxel ( %s, %s, %s ): It appears to have no common-face neighbors "
                                          "inside the image!", str(
                            current_voxel_i), str(current_voxel_j),
                                          str(current_voxel_k))
                        return

                else:
                    out_volume.data[current_voxel_i, current_voxel_j,
                                    current_voxel_k] = labels_nomask[label_index]

        if out_vol_path is None:
            out_vol_path = in_vol_path

        IOUtils.write_volume(out_vol_path, out_volume)

        # Save the output indexes that survived masking
        out_ijk = numpy.vstack(out_ijk)
        filepath = os.path.splitext(out_vol_path)[0]
        numpy.save(filepath + "-idx.npy", out_ijk)
        numpy.savetxt(filepath + "-idx.txt", out_ijk, fmt='%d')

    def vol_val_xyz(self, vol: numpy.ndarray, aff: numpy.ndarray, val: float) -> numpy.ndarray:
        vox_idx = numpy.argwhere(vol == val)
        xyz = aff.dot(numpy.c_[vox_idx, numpy.ones(vox_idx.shape[0])].T)[:3].T
        return xyz

    def compute_label_volume_centers(self, label_volume: numpy.ndarray, affine: numpy.ndarray):

        for val in numpy.unique(label_volume):
            xyz = self.vol_val_xyz(label_volume, affine, val)
            x, y, z = xyz.mean(axis=0)
            yield val, (x, y, z)

    def label_with_dilation(self, to_label_nii_fname: os.PathLike, dilated_nii_fname: os.PathLike,
                            out_nii_fname: os.PathLike):
        """
        Labels a volume using its labeled dilation. The dilated volume is labeled using scipy.ndimage.label function.
        :param to_label_nii_fname: usually a CT-mask.nii.gz
        :param dilated_nii_fname: dilated version of the to_label_nii_fname volume
        """

        # TODO could make dilation with ndimage also.
        mask = IOUtils.read_volume(to_label_nii_fname)
        dil_mask = IOUtils.read_volume(dilated_nii_fname)

        lab, n = scipy.ndimage.label(dil_mask.data)

        # TODO: this change is from tvb-make. Keep it or not? It returns a different result than the old version.
        lab_xyz = list(self.compute_label_volume_centers(lab, dil_mask.affine_matrix))
        lab_sort = numpy.r_[:n + 1]
        # sort labels along AP axis
        for i, (val, _) in enumerate(sorted(lab_xyz, key=lambda t: t[1][1])):
            lab_sort[val] = i
        lab = lab_sort[lab]

        mask.data *= lab
        self.logger.info(
            '%d objects found when labeling the dilated volume.', n)

        IOUtils.write_volume(out_nii_fname, mask)

    def _label_config(self, aparc: nibabel.nifti1.Nifti1Image) -> nibabel.nifti1.Nifti1Image:
        unique_data = numpy.unique(aparc.data)
        unique_data_map = numpy.r_[:unique_data.max() + 1]
        unique_data_map[unique_data] = numpy.r_[:unique_data.size]
        aparc.data = unique_data_map[aparc.data]
        return aparc

    def simple_label_config(self, in_aparc_path: os.PathLike, out_volume_path: os.PathLike):
        """
        Relabel volume to have contiguous values like Mrtrix' labelconfig.
        :param in_aparc_path: volume voxel value is the index of the region it belongs to.
        :return: writes the labeled volume to out_volume_path.
        """

        aparc = IOUtils.read_volume(in_aparc_path)
        aparc = self._label_config(aparc)
        IOUtils.write_volume(out_volume_path, aparc)

    def _label_volume(self, tdi_volume: nibabel.nifti1.Nifti1Image, lo: float=0.5) -> nibabel.nifti1.Nifti1Image:
        mask = tdi_volume.data > lo
        tdi_volume.data[~mask] = 0
        tdi_volume.data[mask] = numpy.r_[1:mask.sum() + 1]
        return tdi_volume

    def label_vol_from_tdi(self, tdi_volume_path: os.PathLike, out_volume_path: os.PathLike, lo: float=0.5):
        """
        Creates a mask of the voxels with tract ends > lo and any other voxels become 0.
        Labels each voxel different from 0 with integer labels starting from 1.
        :param tdi_volume_path: volume voxel value is the sum of tract ends. Voxel without tract ends has value 0.
        :param lo: tract ends threshold used for masking.
        :return: writes labeled volume to :ut_volume_path.
        """

        nii_volume = IOUtils.read_volume(tdi_volume_path)
        tdi_volume = self._label_volume(nii_volume, lo)
        IOUtils.write_volume(out_volume_path, tdi_volume)

    def remove_zero_connectivity_nodes(self, node_volume_path: os.PathLike, connectivity_matrix_path: os.PathLike,
                                       tract_length_path: Optional[str]=None):
        """
        It removes network nodes with zero connectivity from the volume and connectivity matrices.
        The zero connectivity nodes will be labeled with 0 in the volume and the remaining labels will be updated.
        The connectivity matrices will be symmetric.
        :param node_volume_path: tdi_lbl.nii volume path
        :param connectivity_matrix_path: .csv file, output of Mrtrix3 tck2connectome
        :param tract_length_path: optional .csv tract lengths matrix
        :return: overwrites the input volume and matrices with the processed ones. Also saves matrices as .npy.
        """

        node_volume = IOUtils.read_volume(node_volume_path)

        connectivity = numpy.array(numpy.genfromtxt(
            connectivity_matrix_path, dtype='int64'))
        connectivity = connectivity + connectivity.T
        connectivity_row_sum = numpy.sum(connectivity, axis=0)

        nodes_to_keep_indices = connectivity_row_sum > 0
        connectivity = connectivity[nodes_to_keep_indices, :][
                       :, nodes_to_keep_indices]

        numpy.save(os.path.splitext(connectivity_matrix_path)
                   [0] + NPY_EXTENSION, connectivity)
        numpy.savetxt(connectivity_matrix_path, connectivity, fmt='%1d')

        if os.path.exists(str(tract_length_path)):
            connectivity = numpy.array(numpy.genfromtxt(
                tract_length_path, dtype='int64'))
            connectivity = connectivity[nodes_to_keep_indices, :][
                           :, nodes_to_keep_indices]

            numpy.save(os.path.splitext(tract_length_path)
                       [0] + NPY_EXTENSION, connectivity)
            numpy.savetxt(tract_length_path, connectivity, fmt='%1d')

        else:
            self.logger.warning("Path %s is not valid.", tract_length_path)

        nodes_to_remove_indices, = numpy.where(~nodes_to_keep_indices)
        nodes_to_remove_indices += 1

        for node_index in nodes_to_remove_indices:
            node_volume.data[node_volume.data == node_index] = 0

        node_volume.data[node_volume.data > 0] = numpy.r_[
                                                 1:(connectivity.shape[0] + 1)]

        IOUtils.write_volume(node_volume_path, node_volume)

    def con_vox_in_ras(self, ref_vol_path: os.PathLike) -> (numpy.ndarray, numpy.ndarray):
        """
        This function reads a tdi_lbl volume and returns the voxels that correspond to connectome nodes,
        and their coordinates in ras space, simply by applying the affine transform of the volume
        :param ref_vol_path: the path to the tdi_lbl volume
        :return: vox and voxxyz,
                i.e., the labels (integers>=1) and the coordinates of the connnectome nodes-voxels, respectively
        """
        # Read the reference tdi_lbl volume:
        vollbl = IOUtils.read_volume(ref_vol_path)
        vox = vollbl.data.astype('i')
        # Get only the voxels that correspond to connectome nodes:
        voxijk, = numpy.where(vox.flatten() > 0)
        voxijk = numpy.unravel_index(voxijk, vollbl.dimensions)
        vox = vox[voxijk[0], voxijk[1], voxijk[2]]
        # ...and their coordinates in ras xyz space
        voxxzy = vollbl.affine_matrix.dot(numpy.c_[voxijk[0], voxijk[1], voxijk[
            2], numpy.ones(vox.shape[0])].T)[:3].T
        return vox, voxxzy

    def change_labels_of_aparc_aseg(self, atlas_suffix, volume, mapping_dict, conn_regs_nr):
        if atlas_suffix == AtlasSuffix.A2009S:
            volume.data[volume.data == 1000] = 11100
            volume.data[volume.data == 2000] = 12100
        not_matched = set()
        for i in range(volume.data.shape[0]):
            for j in range(volume.data.shape[1]):
                for k in range(volume.data.shape[2]):
                    val = volume.data[i][j][k]
                    if not val in mapping_dict:
                        not_matched.add(val)
                    volume.data[i][j][k] = mapping_dict.get(val, -1)

        print("Now values are in interval [%d - %d]" % (volume.data.min(), volume.data.max()))

        if not_matched:
            print("Not matched regions will be considered background: %s" % not_matched)
        assert (volume.data.min() >= -1 and volume.data.max() < conn_regs_nr)

        return volume

    def transform_coords(self, coords: Union[list, numpy.ndarray, str], src_img: os.PathLike, dest_img: os.PathLike,
                  transform_mat: os.PathLike, output_file: Optional[str]=None) \
            -> (numpy.array, Union[numpy.ndarray, type(None)]):

        if os.path.isfile(coords):
            command = "img2imgcoord %s -mm -src %s -dest %s -xfm %s"
            args = [coords, src_img, dest_img, transform_mat]
        else:
            coords_str = " ".join([str(x) for x in coords])
            command = "echo %s | img2imgcoord -mm -src %s -dest %s -xfm %s"
            args = [coords_str, src_img, dest_img, transform_mat]

        if os.path.isdir(os.path.dirname(output_file)):
            command += " > %s"
            args.append(output_file)

        execute_command(command % tuple(args), cwd=os.path.dirname(transform_mat), shell=True)

        transformed_coords = numpy.loadtxt(output_file, skiprows=1)

        return transformed_coords, output_file
