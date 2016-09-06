#!/bin/bash

# workflow for identifying electrode contact positions from CT and FS recon
# comments with XXX below are visual inspections to do

# script arguments
electrode_intensity_threshold=2000 # magic value, may need adjustment

# check arguments
if [ -z "$SOURCE_CT" ]
then
    echo "SOURCE_CT=/path/to/ct seeg-ct.sh"
    exit 1
fi

# set up topic folder insider subject folder
seeg_dir=$SUBJECTS_DIR/$SUBJECT/seeg
mkdir -p $seeg_dir
pushd $seeg_dir

# convert source CT to standard format (input could be e.g. DICOM folder)
mri_convert $SOURCE_CT CT.nii.gz

# TODO rewrite slightly for resized T1 usage
# TODO set up testing on jenkins

# convert FS images to Nifti for use with FSL
for image in T1 brain aparc+aseg
do
    mri_convert ../mri/$image.mgz $image.nii.gz
done

# standardize image orientations
for image in T1 CT brain aparc+aseg
do
    fslreorient2std $image.nii.gz $image-reo.nii.gz
done

# align to original res first
# resize and applyxfm to CT into resized image
# register CT to T1 space
flirt -in CT-reo.nii.gz \
    -ref  T1-big.nii.gz \
    -omat CT-to-T1.mat \
    -out  CT-in-T1.nii.gz \
    -dof  12 \
    -searchrx -180 180 \
    -searchry -180 180 \
    -searchrz -180 180 \
    -cost mutualinfo \
    -searchcost mutualinfo

# XXX check registration, should fit
freeview -v T1-reo.nii.gz CT-in-T1.nii.gz:opacity=0.5:colormap=jet \
    -viewport sagittal -layout 1 -screenshot ss00-CT-in-T1.png

# invert CT to T1 transform
convert_xfm -omat T1-to-CT.mat -inverse CT-to-T1.mat

# move FS images to CT space (mainly brain mask, others for viz)
for image in T1 brain aparc+aseg
do
    opts="-applyxfm -init T1-to-CT.mat -ref CT-reo.nii.gz"
    # for the volume region mapping, turn off interpolation:
    if [[ $image == "aparc+aseg" ]]
    then
        opts="$opts -interp nearestneighbour"
    fi
    # apply transform
    flirt $opts -in $image-reo.nii.gz -out $image-in-CT.nii.gz
done

# XXX check brain on CT, should fit
freeview -v CT-reo.nii.gz brain-in-CT.nii.gz:opacity=0.5:colormap=jet \
    -viewport axial -layout 1 -screenshot ss01-brain-in-CT.png

# binarize brain with erosion to avoid false positives
mri_binarize --erode 8 \
    --i brain-big.nii \
    --o brain-mask.nii.gz \
    --min 10

# XXX check mask on CT, ovals should line up
freeview -v CT-reo.nii.gz brain-mask.nii.gz:opacity=0.5:colormap=jet \
    -viewport axial -layout 1 -screenshot ss02-brain-mask.png

# apply mask to ct
mri_binarize --i CT-in-T1.nii.gz --o CT-mask.nii.gz \
    --min $electrode_intensity_threshold \
    --mask brain-mask.nii.gz

# check CT electrode mask, if luck otherwise do it interactively to be sure
# XXX CT-mask should highlight the electrodes only, neatly separated
# TODO navigate to axial slice where mask most nonzeros
freeview -v CT-reo.nii.gz CT-mask.nii.gz:opacity=0.5:colormap=jet \
    -viewport axial -layout 1 -screenshot ss03-ct-mask.png

# dilate mask for labeling
mri_binarize --i CT-mask.nii.gz --o CT-dil-mask.nii.gz \
    --min 0.5 --dilate 4 --erode 2

# XXX check dilation (TODO autonav best slice)
freeview -v CT-reo.nii.gz CT-mask.nii.gz:opacity=0.4:colormap=jet \
    CT-dil-mask.nii.gz:colormap=hot:opacity=0.4 \
    -viewport axial -screenshot ss04-ct-dil-mask.png

# label dil mask and apply to mask
# XXX reported number of electrodes should match implantation schema
python<<EOF
import utils
utils.label_with_dilation('CT-mask.nii.gz', 'CT-dil-mask.nii.gz', 'CT-lab-mask.nii.gz')
EOF

# XXX check that labels have different colors, in mrview if possible
# TODO use nibabel to generate subplots per electrode
freeview -v CT-reo.nii.gz CT-lab-mask.nii.gz:opacity=0.5:colormap=jet \
    -viewport axial -screenshot ss05-ct-dil-mask.png

# find contact positions
python<<EOF
import nibabel, numpy as np
nii = nibabel.load('CT-lab-mask.nii.gz')
lab_bin = nii.get_data()
aff = nii.affine
ulab = np.unique(lab_bin)
ulab = ulab[ulab > 0]
ul_to_pos = {}

for ul in ulab[ulab > 0]:
    xyz_pos = periodic_xyz_for_object(lab_bin, ul)
    print(xyz_pos)

    break
    # move to T1 space
    CT_to_T1 = np.loadtxt("CT-to-T1.mat")
    # not correct!? ascii mat is wrong, wtf!
    xyz_T1 = CT_to_T1.dot(np.c_[xyz_pos, np.ones(len(xyz_pos))].T)[:3].T
    # TODO sort these median to lateral
# visual check
if False:
    clf()
    subplot(2, 1, 1)
    plot(bxi, bn)
    subplot(2, 1, 2)
    plot(w, np.abs(Bf))
    subplot(2, 1, 1)
    cos_arg = 2*pi*f[i_peak]*bxi + theta
    plot(bxi, np.cos(cos_arg)*bn.std()+bn.mean(), 'k--', alpha=0.5)
    [axvline(xp, color='r') for xp in xi_pos];
    show()
EOF



# TODO multispaced electrodes like oblique should have weird spectra
# TODO rm electrode artifact in acquisition plane by deconv?

# vim: ts=8 sts=4 sw=4 et ai
