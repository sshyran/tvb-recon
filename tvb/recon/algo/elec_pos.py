import numpy as np
import os
import subprocess
import nibabel as nib


# This is just a helper function I use to run command line stuff
# It is not needed in the workflow that can directly run system commands
def execute_command(command, cwd=os.getcwd(), shell=True):
    import time
    import sys
    print("Running process in directory:\n" + cwd)
    print("Command:\n" + command)
    tic = time.time()
    process = subprocess.Popen(command, shell=shell, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               universal_newlines=True)
    output = []
    while True:
        nextline = process.stdout.readline()
        if nextline == '' and process.poll() is not None:
            break
        sys.stdout.write(nextline)
        sys.stdout.flush()
        output.append(nextline)
    output = "\n".join(output)
    # output = process.communicate()[0]
    exit_code = process.returncode
    if exit_code == 0:
        return output, sys.stdout, time.time() - tic
    else:
        print("The process ran for " + str(time.time() - tic))
        raise subprocess.CalledProcessError(exit_code, command)

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


# I don't know what this is... I have to ask Viktor...
def add_min_max(volume):
    """Ugly hack to deal with the MRtrix bug (?) that causes MRview to crop min/max values"""

    volume[0, 0, 0] = np.nanmin(volume) - 1
    volume[0, 0, 1] = np.nanmax(volume) + 1


# This function creates a new volume matrix from an input label volume matrix
# by setting values[i] to all voxels where input_volume == i
def gen_volume_regions(values, label_volume):
    new_volume = np.zeros(label_volume.shape)
    new_volume[:, :] = np.nan

    for i, value in enumerate(values):
        region = i + 1
        mask = label_volume[:, :, :] == region
        new_volume[mask] = value

    return new_volume


# This function creates a new volume matrix of similar shape to a reference volume volume matrix
# by setting input values at input positions (optionally + & - (int) dist)
# after applying to them the ref_aff transform
def gen_volume_points(values, positions, ref_volume, ref_aff, dist=0):
    new_volume = np.zeros(ref_volume.shape)
    new_volume[:, :] = np.nan

    kx, ky, kz = np.mgrid[-dist:dist + 1, -dist:dist + 1, -dist:dist + 1]

    for val, pos in zip(values, positions):
        ix, iy, iz = np.linalg.solve(ref_aff, np.append(pos, 1.0))[0:3].astype(int)
        # new_volume[inds[0], inds[1], inds[2]] = val
        new_volume[ix + kx, iy + ky, iz + kz] = val

    return new_volume


# Create and save a new nifti label volume with values[i] to all voxels where input_volume == i
def save_nifti_regions(values, label_volume_file, out_file):
    label_nii = nib.load(label_volume_file)
    label_volume = label_nii.get_data()
    new_volume = gen_volume_regions(values, label_volume)
    add_min_max(new_volume)
    new_nii = nib.Nifti1Image(new_volume, label_nii.affine)
    nib.save(new_nii, out_file)


# Create and save a new nifti label volume of similar shape to a reference volume
# by setting input values at input positions (optionally + & - (int) dist)
# after applying to them the ref_aff transform
def save_nifti_points(values, names, position_file, ref_volume_file, out_file, skip_missing=False, dist=0):
    ref_nii = nib.load(ref_volume_file)

    contact_names = list(np.genfromtxt(position_file, dtype=str, usecols=(0,)))
    contact_positions = np.genfromtxt(position_file, dtype=float, usecols=(1, 2, 3))

    positions = np.zeros((len(values), 3))
    for i, name in enumerate(names):
        try:
            contact_ind = contact_names.index(name)
            pos = contact_positions[contact_ind, :]
        except ValueError:
            pos = np.array([np.nan, np.nan, np.nan])

        positions[i, :] = pos

    missing_mask = np.isnan(positions[:, 0])
    if skip_missing:
        values = values[~missing_mask]
        positions = positions[~missing_mask, :]
    else:
        if np.any(missing_mask):
            raise ValueError("Missing contact position(s) for: %s." % ", ".join([
                name for name, missing in zip(names, missing_mask) if missing]))

        # np.array(names)[missing_mask]))

    new_volume = gen_volume_points(values, positions, ref_nii.get_data(), ref_nii.affine, dist=dist)
    add_min_max(new_volume)

    new_nii = nib.Nifti1Image(new_volume, ref_nii.affine)
    nib.save(new_nii, out_file)


# Read the seeg contacts' names and positions from pom file and save then in text and nifti files:
def read_write_pom_files(pomfile, elec_file_pom, mrielec, elec_nii_pom):
    with open(pomfile) as fl:
        lines = fl.readlines()

    loc_start = find_line_starting_with(lines, "LOCATION_LIST START_LIST") + 1
    loc_end = find_line_starting_with(lines, "LOCATION_LIST END_LIST")

    names_start = find_line_starting_with(lines, "REMARK_LIST START_LIST") + 1
    names_end = find_line_starting_with(lines, "REMARK_LIST END_LIST")

    coords_list = lines[loc_start:loc_end]
    names = lines[names_start:names_end]

    coords_list = np.array([list(map(float, coords.split("\t"))) for coords in coords_list])
    names = [name.strip() for name in names]

    save_xyz_file(coords_list, names, elec_file_pom)
    save_nifti_points(np.ones(len(names)), names, elec_file_pom, mrielec, elec_nii_pom, dist=1)

    return names, coords_list


# TODO: an automatic way to pick contact's coordinates randomly or one per electrode from the seeg_pom file
coords_from = np.array([
    [-20.59, 19.51, 38.49],
    [-71.79, 10.24, 38.51],
    [2.95, 43.21, 105.8],
    [-48.5, 22.58, 104.8],
    [-0.69, 74.26, 26.77],
    [-60.11, 70.55, 30.76],
    [3.822, -10.42, 100.8]
])

# TODO: an automatic way to pick the corresponding contacts'coordinates from the mri_electrode volume
coords_to = np.array([
    [-21.25, 16.31, 29.61],
    [-72.75, 7.59, 28.61],
    [1.882, 29.46, 99.93],
    [-50.33, 9.83, 95.7],
    [-2.35, 71.61, 27.54],
    [-60.82, 68.77, 30.01],
    [1.705, -22.75, 87.75]
])

# Find the affine transformation pomfile -> MRIelectrodes
def pom_in_mrielec_transform(coords_from, coords_to):

    # Pad the data with ones, so that our transformation can do translations too
    n = coords_from.shape[0]
    pad = lambda x: np.hstack([x, np.ones((x.shape[0], 1))])
    unpad = lambda x: x[:, :-1]
    X = pad(coords_from)
    Y = pad(coords_to)

    # Solve the least squares problem X * A = Y
    # to find our transformation matrix A
    A, res, rank, s = np.linalg.lstsq(X, Y)

    aff_transform = lambda x: unpad(np.dot(pad(x), A))

    print("Target:")
    print(coords_to)
    print("Result:")
    print(aff_transform(coords_from))
    print("Max error:", np.abs(coords_to - aff_transform(coords_from)).max())

    return aff_transform


# Transform the coordinates along (pom file ->) mri_elec -> t1
def transform(elec_file_pom, src_img, dest_img, elec_file_t1, elec_nii_t1, transform_mat, aff_transform=None):
    names = []
    coords = []
    with open(elec_file_pom) as fl:
        lines = fl.readlines()
    for line in lines:
        split_line = line.split()
        names.append(split_line[0])
        coords.append([float(s) for s in [split_line[1], split_line[2], split_line[3]]])
    coords = np.array(coords)
    if aff_transform is not None:
        coords = aff_transform(coords.reshape((1, 3))).reshape((3,))
    dirname = os.path.dirname(elec_file_pom)
    coords_file = os.path.join(dirname, "seeg_coords.xyz")
    np.savetxt(coords_file, coords, fmt='%.3e', delimiter=' ')
    coords_file_t1 = os.path.join(dirname, "seeg_coords_in_t1.xyz")
    output, std_out, time = \
        execute_command("img2imgcoord %s -mm -src %s -dest %s -xfm %s > %s" \
                            % (coords_file, src_img, dest_img, transform_mat, coords_file_t1),
                        cwd=os.path.dirname(transform_mat), shell=True)
    # Save final coordinates in t1 space to text and nifti files
    with open(coords_file_t1) as fl:
        lines = fl.readlines()
    newlines = []
    for line, name in zip(lines[1:], names):
        newlines.append(name + " " + line)
    with open(elec_file_t1, "w") as fl:
        fl.writelines(newlines)
    save_nifti_points(np.ones(len(names)), names, elec_file_t1, dest_img, elec_nii_t1, dist=1)
    # # Originally in Viktor's code (but I was getting an error that there is no run attribute for module...:
    # cp = subprocess.run("echo %s | img2imgcoord -mm -src %s -dest %s -xfm %s" \
    #                     % (coords_str, src_img, dest_img, transform_mat),
    #                     shell=True, stdout=subprocess.PIPE)
    # transformed_coords_str = cp.stdout.decode('ascii').strip().split('\n')[-1]
    # transformed_coords =  np.array([float(x) for x in transformed_coords_str.split(" ") if x])



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
    mrielec_to_t1 = os.path.join(elec_folder, "mrielec_to_t1.mat")
    mrielec_in_t1 = os.path.join(elec_folder, "mrielec_in_t1.nii.gz")
    elec_file_t1 = os.path.join(elec_folder, "seeg.xyz")
    elec_nii_t1 = os.path.join(elec_folder, "elec_t1.nii.gz")


    # This is a 2 to 4 step/job process:

    # A. Preprocessing job. Not necessary if something like seeg_pom.xyz was alraedy given in the input
    # Read the seeg contacts' names and positions from pom file and save then in text and nifti files:
    names, coords_list = read_write_pom_files(pomfile, elec_file_pom, mrielec, elec_nii_pom)

    # B. Optional & the ONLY MANUAL job in case the pom coordinates are not in the same space as the MRIelectrode volume
    # Find the affine transformation pomfile -> MRIelectrodes
    if POM_TO_MRIELEC_TRNSFRM:
        aff_transform = pom_in_mrielec_transform()
    else:
        aff_transform = None

    # C. Optional job in case we haven't computed already the MRIelec_in_t1 transformation matrix
    #    (This must be already part of the workflow...)
    if not os.path.isfile(mrielec_to_t1):
        # Align mrielec to t1
        # Options for flirt co-registration
        regopts = "-cost mutualinfo -dof 12 -searchrz -180 180 -searchry -180 180  -searchrx -180 180"
        # regopts = "-searchrz -20 20 -searchry -20 20 -searchrx -20 20"
        output, std_out, time = \
            execute_command("flirt "+ regopts +
                            " -in " + mrielec +
                            " -ref " + t1 +
                            " -omat " + mrielec_to_t1 +
                            " -out " + mrielec_in_t1, cwd=results_root, shell=True)

    # D. Make the actual transformation from pom/mri_elec space to t1 space
    # Transform the coordinates along (pom file ->) mri_elec -> t1 and save the to text and nifti files
    transform(elec_file_pom, mrielec,  t1,   elec_file_t1, elec_nii_t1, mrielec_to_t1, aff_transform)

if __name__ == "__main__":

    main_elec_pos("TVB1", False)
