#!/usr/bin/env bash

source //anaconda/bin/activate tvb_recon_python3_env

export FREESURFER_HOME
export SUBJECTS_DIR
export SUBJECT=$4

labeled_CT="'$1'"
schema="'$2'"
output="'$3'"

python <<EOF
import nibabel, numpy as np, tvb.recon.algo.reconutils as utils
label_to_name = {}
with open($schema, "r") as fd:
    for line in fd.readlines():
        label_num, label_name = line.strip().split(' ')
        label_num = int(label_num)
        assert label_num not in label_to_name
        label_to_name[int(label_num)] = label_name
nii = nibabel.load($labeled_CT)
lab_bin = nii.get_data()
aff = nii.affine
ulab = np.unique(lab_bin)
ulab = ulab[ulab > 0]
# TODO find closest voxel in parcellation, provide name
fmt = '%s%d\t%f\t%f\t%f\n'
with open($output, "w") as fd:
    for ul in ulab[ulab > 0]:
        if ul not in label_to_name:
            print("skipping object label %d" % (ul, ))
            continue
        xyz_pos = utils.periodic_xyz_for_object(lab_bin, ul, aff)
        xyz_pos = xyz_pos[np.argsort(np.abs(xyz_pos[:, 0]))]
        name = label_to_name[ul]
        for i, (x, y, z) in enumerate(xyz_pos):
            fd.write(fmt % (name, i, x, y, z))
EOF