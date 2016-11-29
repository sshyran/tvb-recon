#!/bin/bash

#Post recon-all:

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



