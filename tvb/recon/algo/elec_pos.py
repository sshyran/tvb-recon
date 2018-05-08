import numpy as np
import os
import nibabel
from tvb.recon.algo.service.utils import execute_command  # , compute_affine_transform
from tvb.recon.algo.service.volume import VolumeService
from tvb.recon.algo.service.sensor import SensorService


volume_service = VolumeService()
sensor_service = SensorService()


# Search a line in a file starting with a string pattern
def find_line_starting_with(line_list, string):
    for i, line in enumerate(line_list):
        if line[:len(string)] == string:
            return i
    return None


# Save coordinates to a file
def save_xyz_file(coords, names, filename):
    with open(filename, "w") as fl:
        for name, coord in zip(names, coords):
            fl.write("%s %8.2f %8.2f %8.2f\n" % (name, coord[0], coord[1], coord[2]))


# Read the seeg contacts' labels and positions from pom file and save then in text and nifti files:
def read_write_pom_files(pomfile, elec_file_pom, mrielec, elec_nii_pom):
    """
    # NEEDED ONLY FOR CC PATIENTS!
    # Read the seeg contacts' labels and positions from pom file and save then in text and nifti files:

    :param pomfile:
    :param elec_file_pom:
    :param mrielec:
    :param elec_nii_pom:
    :return:
    """
    with open(pomfile) as fl:
        lines = fl.readlines()

    loc_start = find_line_starting_with(lines, "LOCATION_LIST START_LIST") + 1
    loc_end = find_line_starting_with(lines, "LOCATION_LIST END_LIST")

    labels_start = find_line_starting_with(lines, "REMARK_LIST START_LIST") + 1
    labels_end = find_line_starting_with(lines, "REMARK_LIST END_LIST")

    coords_list = lines[loc_start:loc_end]
    labels = lines[labels_start:labels_end]

    coords_list = np.array([list(map(float, coords.split("\t"))) for coords in coords_list])
    labels = [label.strip() for label in labels]

    save_xyz_file(coords_list, labels, elec_file_pom)
    n_coords = coords_list.shape[0]
    volume_service.gen_label_volume_from_coords(coords_list, mrielec, elec_nii_pom,
                                                values=np.array(range(n_coords)) + 1, labels=labels,
                                                skip_missing=False, dist=1)

    return labels, coords_list


def binarize_dil_erod(in_vol, out_vol, th=1, dilate=0, erode=0):
    command = "mri_binarize " + " --i " + in_vol + " --o " + out_vol + " --min " + str(th)
    if dilate > 0:
        command += " --dilate " + str(dilate)
    if erode > 0:
        command += " --erode " + str(erode)
    execute_command(command)


def extract_seeg_contacts_from_mrielec(mrielec, mrielec_seeg, dilate=5, erode=1):
    volume = nibabel.load(mrielec)
    data = volume.get_data()
    max_voxel = data.max()
    del data
    binarize_dil_erod(mrielec, mrielec_seeg, max_voxel, dilate, erode)


def coregister_elec_pom_and_mri(seeg_pom_vol, mrielec, seeg_pom_xyz, elec_folder=None, dilate=0, erode=0):

    # Binarize mrielec to get only the seeg contacts
    mrielec_seeg = os.path.join(elec_folder, "mrielec_seeg.nii.gz")
    if not os.path.isfile(mrielec_seeg):
        extract_seeg_contacts_from_mrielec(mrielec, mrielec_seeg, dilate, erode)

    # Binarize the pom seeg contacts
    seeg_nii_pom_bin = seeg_pom_vol.split(".")[0] + "_bin.nii.gz"
    if not os.path.isfile(seeg_nii_pom_bin):
        binarize_dil_erod(seeg_pom_vol, seeg_nii_pom_bin, 1, dilate, erode)

    if not(os.path.isdir(elec_folder)):
        elec_folder = os.path.abspath(seeg_pom_vol)

    # Co-register seeg pom and mrielec contacts
    pom_to_mrielec = os.path.join(elec_folder, "pom_to_mrielec.mat")
    pom_in_mrielec_seeg = os.path.join(elec_folder, "pom_in_mrielec_seeg.nii.gz")
    if not os.path.isfile(pom_to_mrielec):
        # This takes time. It should be independent step and check before repeating
        regopts = "-dof 12 -searchrz -180 180 -searchry -180 180  -searchrx -180 180"
        execute_command("flirt " + regopts +
                        " -in " + seeg_nii_pom_bin +
                        # " -inweight " + elec_nii_pom_dil +
                        " -ref " + mrielec_seeg +
                        # " -refweight " + mrielec_dil +
                        " -omat " + pom_to_mrielec +
                        " -out " + pom_in_mrielec_seeg)

    # Generate the new seeg labels and coords file in the mrielec space, as well as a volume for visual checking
    seeg_pom_xyz_in_mrielec = os.path.join(elec_folder, "seeg_pom_xyz_in_mrielec.xyz")
    pom_vol_in_mrielec = os.path.join(elec_folder, "pom_in_mrielec_vol.nii.gz")
    if not os.path.isfile(seeg_pom_xyz_in_mrielec):
        transform(seeg_pom_xyz, seeg_pom_vol, mrielec, seeg_pom_xyz_in_mrielec, pom_vol_in_mrielec, pom_to_mrielec)

    return seeg_pom_xyz_in_mrielec, pom_vol_in_mrielec


# Transform the coordinates and create output contact files (labels + coords) and a volume for visual checking
def transform(seeg_xyz_in, input_vol, ref_vol, seeg_xyz_ref, seeg_vol_ref, transform_mat, aff_transform=None):
    labels, coords = sensor_service.read_seeg_labels_coords_file(seeg_xyz_in)
    if aff_transform is not None:
        coords = aff_transform(coords.reshape((1, 3))).reshape((3,))
    coords_xyz_in = seeg_xyz_in.replace("seeg", "coords")
    np.savetxt(coords_xyz_in, coords, fmt='%.3e', delimiter=' ')
    coords_xyz_ref = seeg_xyz_ref.replace("seeg", "coords")
    coords_in_ref, coords_xyz_ref = volume_service.transform_coords(coords_xyz_in, input_vol, ref_vol, transform_mat,
                                                                    coords_xyz_ref)
    # Save final coordinates in ref_vol space to text and nifti files
    with open(coords_xyz_ref) as fl:
        lines = fl.readlines()
    newlines = []
    for line, label in zip(lines[1:], labels):
        newlines.append(label + " " + line)
    with open(seeg_xyz_ref, "w") as fl:
        fl.writelines(newlines)
    volume_service.gen_label_volume_from_coords(coords_in_ref, ref_vol, seeg_vol_ref, labels=labels,
                                                skip_missing=False, dist=1)


def main_elec_pos(patient, POM_TO_MRIELEC_TRNSFRM=False, dilate=0, erode=0):

    # Paths:

    root = "/Users/dionperd/Dropbox/Work/VBtech/VEP"
    results_root = os.path.join(root, "results/CC")
    data_root = os.path.join(root, "data/CC")
    # patient = "TVB3"

    #Inputs:
    t1 = os.path.join(results_root, patient, "tvb", "T1.nii.gz")
    mrielec = os.path.join(data_root, patient, "raw/mri", "MRIelectrode.nii.gz")
    pomfile = os.path.join(data_root, patient, "raw", "Electrodes.pom")

    # Generated files
    elec_folder = os.path.join(data_root, patient, "elecs")
    if not os.path.isdir(elec_folder):
        os.mkdir(elec_folder)
    seeg_pom_xyz = os.path.join(elec_folder, "seeg_pom_xyz.txt")
    seeg_pom_vol = os.path.join(elec_folder, "seeg_pom_vol.nii.gz")
    mrielec_to_t1 = os.path.join(elec_folder, "mrielec_to_t1.mat") # ELEC_to_T1 or mrielec_to_t1
    mrielec_in_t1 = os.path.join(elec_folder, "mrielec_in_t1.nii.gz") # ELEC_in_T1 or mrielec_in_t1
    seeg_xyz = os.path.join(elec_folder, "seeg_xyz.txt")
    seeg_in_t1 = os.path.join(elec_folder, "seeg_in_t1.nii.gz")


    # This is a 2 to 4 step/job process:

    # A. Preprocessing job. Not necessary if something like seeg_pom.xyz was alraedy given in the input
    # Read the seeg contacts' names and positions from pom file and save then in text and nifti volume files:
    if not os.path.isfile(seeg_pom_vol):
        names, coords_list = read_write_pom_files(pomfile, seeg_pom_xyz, mrielec, seeg_pom_vol)

    # B. Optional job in case we haven't computed already the MRIelec_in_t1 transformation matrix
    #    (This must be already part of the workflow...)
    # This takes time. It should be independent step and check before repeating
    if not os.path.isfile(mrielec_to_t1):
        # Align mrielec to t1
        # Options for flirt co-registration
        regopts = "-cost mutualinfo -dof 12 -searchrz -180 180 -searchry -180 180  -searchrx -180 180"
        # regopts = "-searchrz -20 20 -searchry -20 20 -searchrx -20 20"
        execute_command("flirt " + regopts +
                        " -in " + mrielec +
                        " -ref " + t1 +
                        " -omat " + mrielec_to_t1 +
                        " -out " + mrielec_in_t1)

    # C. Optional job in case the pom coordinates are not in the same space as the MRIelectrode volume
    # Find the affine transformation pomfile -> MRIelectrodes
    if POM_TO_MRIELEC_TRNSFRM:
        # Alternative manual way:
        # This is how Viktor Sip did it. It requires manually specifying the point coordinates below.
        # coords_from = np.array([
        #     [-20.59, 19.51, 38.49],
        #     [-71.79, 10.24, 38.51],
        #     [2.95, 43.21, 105.8],
        #     [-48.5, 22.58, 104.8],
        #     [-0.69, 74.26, 26.77],
        #     [-60.11, 70.55, 30.76],
        #     [3.822, -10.42, 100.8]
        # ])
        #
        # coords_to = np.array([
        #     [-21.25, 16.31, 29.61],
        #     [-72.75, 7.59, 28.61],
        #     [1.882, 29.46, 99.93],
        #     [-50.33, 9.83, 95.7],
        #     [-2.35, 71.61, 27.54],
        #     [-60.82, 68.77, 30.01],
        #     [1.705, -22.75, 87.75]
        # ])
        # aff_transform = compute_affine_transform(coords_from, coords_to)

        seeg_pom_xyz_in_mrielec, pom_vol_in_mrielec = \
            coregister_elec_pom_and_mri(seeg_pom_vol, mrielec,  seeg_pom_xyz, elec_folder, dilate, erode)

        # D. Make the actual transformation from pom/mri_elec space to t1 space
        transform(seeg_pom_xyz_in_mrielec, pom_vol_in_mrielec, t1, seeg_xyz, seeg_in_t1, mrielec_to_t1)

    else:

        # D. Make the actual transformation from pom/mri_elec space to t1 space
        # Transform the coordinates along (pom file ->) mri_elec -> t1 and save the result to text and nifti files
        transform(seeg_pom_xyz, mrielec, t1, seeg_xyz, seeg_in_t1, mrielec_to_t1)


if __name__ == "__main__":

    main_elec_pos("TVB1", True, 0, 0)
