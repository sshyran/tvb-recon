#!/bin/bash

# workflow for identifying electrode contact positions from CT and FS recon
# comments with XXX below are visual inspections to do

if [ -z "$SUBJECT" || -z "$CT" ]
then
    echo "CT=/path/to/ct/data SUBJECT=fs_subject_name seeg-ct.sh"
    exit 1
fi

set -eu
set -o pipefail

threshold=${CT_ELEC_INTENSITY_TH:-"2000"}

topic_folder=$SUBJECTS_DIR/$SUBJECT/seeg
mkdir -p $topic_folder
pushd $topic_folder

# input CT to nii gz
mri_convert $CT ./CT.nii.gz --out_orientation RAS -rt nearest

# fslreorient2std ./CT.nii.gz ./CT-reo.nii.gz

for image in T1 brain ; do
    # fslreorient2std $MRI/$image.nii.gz ./$image-reo.nii.gz
    mri_convert ../mri/$image.mgz ./$image-big.nii.gz -ds 0.5 0.5 0.5
done

mri_convert ../mri/T1.mgz ./T1.nii.gz --out_orientation RAS
flirt -v -in ./CT.nii.gz -ref ./T1.nii.gz -omat CT_to_T1.mat -out CT_in_T1.nii.gz \
    -dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -cost mutualinfo

# binarize brain mask; erode to avoid inner skull
mri_binarize --i ../mri/brain.mgz --o mask.nii --min 10 --erode 8

# mask ct
mri_binarize --i CT_in_T1.nii.gz --o masked_CT.nii --min $threshold --mask mask.nii

# dilate mask for labeling
mri_binarize --i masked_CT.nii --o dilated_CT.nii --min 0.5 --dilate 4 --erode 2

# label dilated mask and apply to mask
python3 <<EOF
from tvb.recon.algo.reconutils import label_with_dilation
label_with_dilation("masked_CT.nii", "dilated_CT.nii", "labeled_CT.nii")
EOF

# TODO compare labels to implantation schema
mrview T1.nii.gz -overlay.load=dilated_CT.nii
$open $implantation_pptx

# TODO and put correspondence in file
cat > schema.txt <<EOF
0 FCA'
1 GCA'
etc
EOF
# labels w/o names will be ignored

# find contact positions
python<<EOF
import nibabel, numpy as np, tvb.recon.algo.reconutils as utils
label_to_name = {}
with open("schema.txt", "r") as fd:
    for line in fd.readlines():
        label_num, label_name = line.strip().split(' ')
        label_to_name[int(label_num)] = label_name
nii = nibabel.load("")
lab_bin = nii.get_data()
aff = nii.affine
ulab = np.unique(lab_bin)
ulab = ulab[ulab > 0]
# TODO find closest voxel in parcellation, provide name
fmt = '%s%d\t%f\t%f\t%f\n'
with open("$seeg_positions", "w") as fd:
    for ul in ulab[ulab > 0]:
        xyz_pos = utils.periodic_xyz_for_object(lab_bin, ul, aff)
        name = label_to_name[ul]
        for i, (x, y, z) in enumerate(xyz_pos):
            fd.write(fmt % (i, x, y, z))
EOF

