#!/bin/bash

#Post recon-all:

pushd $MRI

#Get CRAS into an environment variable and into a text file
CRAS="$(mri_info --cras ./T1.mgz)"
mri_info --cras ./T1.mgz > $CRAS_PATH

#Get vox2ras and vox2tkras into text files
mri_info --vox2ras ./T1.mgz > $T1_NAT_VOX2RAS_PATH
mri_info --vox2ras-tkr ./T1.mgz > $T1_NAT_VOX2RASTKR_PATH
mri_info --vox2ras ./T1.nii.gz > $T1_VOX2RAS_PATH
mri_info --vox2ras-tkr ./T1.nii.gz > $T1_VOX2RASTKR_PATH

#Generate nifti files with good orientation
for vol in aparc+aseg aseg norm orig wm
do
    mri_convert $MRI/$vol.mgz $MRI/$vol.nii.gz --out_orientation RAS -rt nearest
    #!!Probably not necessary any more!:
    #fslreorient2std $MRI/$vol.nii.gz $MRI/$vol-reo.nii.gz
    #mv $MRI/$vol-reo.nii.gz $MRI/$vol.nii.gz
done

#!!Probably not necessary any more!:
#Get the wm surface for the cortical regions:
#for h in lh rh
#do
#    #Transform from surface tkRAS coordinates to scanner RAS ones:
#    mris_convert --to-scanner $SURF/$h.white $SURF/$h.white-ras
#    #Create a white surface annotation
#    cp $LABEL/$h.$DEFAULT_APARC.annot $LABEL/$h.white.annot
#done



