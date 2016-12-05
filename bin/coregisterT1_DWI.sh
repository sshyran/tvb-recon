#!/bin/bash

#FLIRT co-registration of T1 with DWI

pushd $DMR

if [ "$COREG_USE" = "flirt" ]
then
    #Register DWI to T1 and get the relevant transform
    regopt="-dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -cost mutualinfo"
    flirt -in ./b0.nii.gz -ref $MRI/T1.nii.gz -omat ./d2t.mat -out ./b0-in-t1.nii.gz $regopt

    #Generate and apply the inverse transform from T1 to DWI for T1
    convert_xfm -omat ./t2d.mat -inverse ./d2t.mat
    flirt -applyxfm -in $MRI/T1.nii.gz -ref ./b0.nii.gz -out ./t1-in-d.nii.gz -init ./t2d.mat
else
    #Register DWI to T1 and get the relevant transform
    bbregister --s $SUBJECT --mov ./b0.nii.gz --o ./b0-in-t1.mgz --dti --reg ./d2t.reg --lta ./d2t.lta --fslmat ./d2t.mat

    #Convert to ras coordinates
    mri_convert ./b0-in-t1.mgz ./b0-in-t1.nii.gz --out_orientation ras
    rm ./b0-in-t1.mgz

    #Apply the inverse transform from T1 to DWI for T1
    mri_vol2vol --mov $MRI/T1.mgz --targ ./b0.nii.gz --o ./t1-in-d.nii.gz --lta-inv ./d2t.lta --save-reg
    mv ./t1-in-d.nii.gz.lta ./t2d.lta
    mv ./t1-in-d.nii.gz.reg ./t2d.reg
fi

#Visual check:
#mrview ./b0.nii.gz -overlay.load ./T1-in-d.nii.gz -overlay.opacity 0.3 -mode 2 #-mode 2 is the view
#source $CODE/snapshot.sh use_freeview 2vols ./b0.nii.gz ./t1-in-d.nii.gz
python -m $SNAPSHOT --snapshot_name b0_t1_in_d --ras_transform 2vols ./b0.nii.gz ./t1-in-d.nii.gz
python -m $SNAPSHOT --snapshot_name b0_t1_in_t1 2vols $MRI/T1.nii.gz ./b0-in-t1.nii.gz

popd




