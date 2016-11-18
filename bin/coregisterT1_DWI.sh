#!/bin/bash

#FLIRT co-registration of T1 with DWI

pushd $DMR

#Convert T1 to NIFTI with good orientation
mri_convert $MRI/T1.mgz $MRI/T1.nii.gz --out_orientation RAS -rt nearest
fslreorient2std $MRI/T1.nii.gz $MRI/T1-reo.nii.gz
mv $MRI/T1-reo.nii.gz $MRI/T1.nii.gz

#Register DWI to T1 and get the relevant transform
regopt="-dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -cost mutualinfo"
flirt -in ./b0.nii.gz -ref $MRI/T1.nii.gz -omat ./d2t.mat -out ./b0-in-t1.nii.gz $regopt

#Generate and apply the inverse transform from T1 to DWI for T1
convert_xfm -omat ./t2d.mat -inverse ./d2t.mat
flirt -applyxfm -in $MRI/T1.nii.gz -ref ./b0.nii.gz -out ./T1-in-d.nii.gz -init ./t2d.mat -interp nearestneighbour #the nn algorithm is needed because we have integer values

#Visual check:
#mrview ./b0.nii.gz -overlay.load ./T1-in-d.nii.gz -overlay.opacity 0.3 -mode 2 #-mode 2 is the view
source $CODE/snapshot.sh use_freeview 2vols ./b0.nii.gz ./T1-in-d.nii.gz

popd




