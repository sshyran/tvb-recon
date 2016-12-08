#!/bin/bash

#FLIRT co-registration of gmwmi with aseg

pushd $ASEG_VOL

#Anatomically constraint spherical deconvolution
mrconvert $DMR/tdi_ends.mif ./tdi_ends.nii

mri_convert $MRI/aseg.mgz $MRI/aseg.nii.gz --out_orientation RAS -rt nearest
fslreorient2std $MRI/aseg.nii.gz $MRI/aseg-reo.nii.gz
mv $MRI/aseg-reo.nii.gz $MRI/aseg.nii.gz

flirt -applyxfm -in $MRI/aseg.nii.gz -ref $DMR/b0.nii.gz -out ./aseg-in-d.nii.gz -init ./t2d.mat -interp nearestneighbour


#Scale gmwmi-in-T1 in 0-256 for visualization reasons
mris_calc ./gmwmi-in-T1.nii.gz mul 256
mri_convert ./out.mgz ./gmwmi-in-T1-256.nii.gz
rm ./gmwmi-in-T1.nii.gz

#Convert aseg to NIFTI with good orientation

#Visual checks
#(interactive):
freeview -v $MRI/T1.nii.gz ./aseg.nii.gz ./gmwmi-in-T1-256.nii.gz:opacity=0.5
#(screenshot):
freeview -v $MRI/T1.nii.gz ./aseg.nii.gz ./gmwmi-in-T1-256.nii.gz:opacity=0.5 -ss $FIGS/gmwmi-in-T1-aseg.png


#Construct surfaces around the original aseg and mask them with gmwmi for each label separately
for l in $ASEG_LIST
do
    lbl=$(printf %03d $l)

    #Construct surface
    sh $CODE/aseg2srf_DP.sh -s $SUBJECT -l $l -d

    for vn in 0 1 2 3 4 5
    do
        #Mask with gmwmi and save the resulting surface in a new file
        python -c "import bnm.recon.algo.reconutils; bnm.recon.algo.reconutils.gmwmi_to_surf('./aseg_$lbl','./gmwmi-in-T1.nii.gz','./aseg_$lbl-gw-$vn',vn=$vn,th=$ASEG_GW_THR)"

    done
done
#Visual check (screenshot):
freeview -f ./aseg*:color=grey ./gmwmi-in-T1.nii.gz:opacity=0.5 -viewport 3d -ss $FIGS/gmwmi-aseg-surfs.png

popd



