import numpy as np
import os
from tvb.recon.algo.service.utils import compute_affine_transform, execute_command
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
    volume_service.gen_label_volume_from_coords(np.ones(len(labels)), np.array(coords_list), labels, mrielec,
                                                elec_nii_pom, skip_missing=False, dist=1)

    return labels, coords_list


# Transform the coordinates along (pom file ->) mri_elec -> t1
def transform(elec_file_pom, mri_elec, t1, elec_file_t1, elec_nii_t1, transform_mat, aff_transform=None):
    labels, coords = sensor_service.read_seeg_labels_coords_file(elec_file_pom)
    if aff_transform is not None:
        coords = aff_transform(coords.reshape((1, 3))).reshape((3,))
    dirname = os.path.dirname(elec_file_pom)
    coords_file = os.path.join(dirname, "seeg_coords.xyz")
    np.savetxt(coords_file, coords, fmt='%.3e', delimiter=' ')
    coords_file_t1 = os.path.join(dirname, "seeg_coords_in_t1.xyz")
    # TODO correct it to get coords_in_t1
    coords_in_t1, coords_file_t1 = volume_service.transform_coords(coords_file, mri_elec, t1, transform_mat,
                                                                   coords_file_t1)
    # output, std_out, time = \
    #     execute_command("img2imgcoord %s -mm -src %s -dest %s -xfm %s > %s" \
    #                         % (coords_file, src_img, dest_img, transform_mat, coords_file_t1),
    #                     cwd=os.path.dirname(transform_mat), shell=True)
    # Save final coordinates in t1 space to text and nifti files
    with open(coords_file_t1) as fl:
        lines = fl.readlines()
    newlines = []
    for line, label in zip(lines[1:], labels):
        newlines.append(label + " " + line)
    with open(elec_file_t1, "w") as fl:
        fl.writelines(newlines)
    volume_service.gen_label_volume_from_coords(np.ones(len(labels)), elec_file_t1, labels, t1, elec_nii_t1,
                                                skip_missing=False, dist=1)
    # volume_service.gen_label_volume_from_coords(np.ones(len(labels)), coords_in_t1, labels, t1, elec_nii_t1,
    #                                             skip_missing=False, dist=1)


def main_elec_pos(patient, POM_TO_MRIELEC_TRNSFRM=False):

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
    elec_file_pom = os.path.join(elec_folder, "seeg_pom.xyz")
    elec_nii_pom = os.path.join(elec_folder, "elec_pom.nii.gz")
    mrielec_to_t1 = os.path.join(elec_folder, "ELEC_to_T1.mat") # ELEC_to_T1 or mrielec_to_t1
    mrielec_in_t1 = os.path.join(elec_folder, "ELEC_in_T1.nii.gz") # ELEC_in_T1 or mrielec_in_t1
    elec_file_t1 = os.path.join(elec_folder, "seeg.xyz")
    elec_nii_t1 = os.path.join(elec_folder, "elec_t1.nii.gz")


    # This is a 2 to 4 step/job process:

    # A. Preprocessing job. Not necessary if something like seeg_pom.xyz was alraedy given in the input
    # Read the seeg contacts' names and positions from pom file and save then in text and nifti files:
    names, coords_list = read_write_pom_files(pomfile, elec_file_pom, mrielec, elec_nii_pom)

    # B. Optional & the ONLY MANUAL job in case the pom coordinates are not in the same space as the MRIelectrode volume
    # Find the affine transformation pomfile -> MRIelectrodes
    if POM_TO_MRIELEC_TRNSFRM:
        # # TODO: an automatic way to pick contact's coordinates randomly or one per electrode from the seeg_pom file
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
        # # TODO: an automatic way to pick the corresponding contacts'coordinates from the mri_electrode volume
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

        # Alternative way:
        pom_to_mrielec = os.path.join(elec_folder, "pom_to_mrielec.mat")
        pom_in_mrielec = os.path.join(elec_folder, "pom_in_mrielec.nii.gz")

        regopts = "-cost mutualinfo -dof 12 -searchrz -180 180 -searchry -180 180  -searchrx -180 180"
        execute_command("flirt " + regopts +
                        " -in " + elec_nii_pom +
                        " -inweight " + elec_nii_pom +
                        " -ref " + mrielec +
                        " -refweight " + mrielec +
                        " -omat " + pom_to_mrielec +
                        " -out " + pom_in_mrielec)

        elec_file_mrielec = os.path.join(elec_folder, "seeg_pom_in_mrielec.xyz")
        transform(elec_file_pom, elec_nii_pom, mrielec, elec_file_mrielec, pom_in_mrielec, pom_to_mrielec)

    else:
        aff_transform = None

    #
    # C. Optional job in case we haven't computed already the MRIelec_in_t1 transformation matrix
    #    (This must be already part of the workflow...)
    if not os.path.isfile(mrielec_to_t1):
        # Align mrielec to t1
        # Options for flirt co-registration
        regopts = "-dof 12 -searchrz -180 180 -searchry -180 180  -searchrx -180 180"
        # regopts = "-searchrz -20 20 -searchry -20 20 -searchrx -20 20"
        execute_command("flirt " + regopts +
                        " -in " + mrielec +
                        " -ref " + t1 +
                        " -omat " + mrielec_to_t1 +
                        " -out " + mrielec_in_t1)

    # D. Make the actual transformation from pom/mri_elec space to t1 space
    # Transform the coordinates along (pom file ->) mri_elec -> t1 and save the result to text and nifti files
    transform(elec_file_pom, mrielec, t1, elec_file_t1, elec_nii_t1, mrielec_to_t1, aff_transform)

    # # Alternative way:
    #
    # execute_command("flirt " + regopts +
    #                 " -in " + pom_in_mrielec +
    #                 " -inweight " + pom_in_mrielec +
    #                 " -ref " + t1 +
    #                 " -omat " + mrielec_to_t1 +
    #                 " -out " + elec_nii_t1)
    # transform(elec_file_mrielec, pom_in_mrielec, t1, elec_file_t1, elec_nii_t1, mrielec_to_t1)

if __name__ == "__main__":

    main_elec_pos("TVB4", True)
