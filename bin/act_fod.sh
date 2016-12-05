#!/bin/bash

pushd $DMR

#Anatomically constraint spherical deconvolution
5ttgen fsl ./t1-in-d.nii.gz ./5tt.mif  #if a brain mask is already applied: -premasked
5tt2gmwmi ./5tt.mif ./gmwmi.mif -nthreads $MRTRIX_THRDS -force

#Visual checks (interactive):
#Not really worthwhile:
5tt2vis ./5tt.mif ./5ttvis.mif -force
mrview ./5ttvis.mif
#More useful:
#-mode 2 is the view
mrview ./b0.nii.gz -overlay.load ./5ttvis.mif -overlay.opacity 0.3 -mode 2 &
mrconvert ./5ttvis.mif ./5ttvis.nii.gz -force
#source snapshot.sh use_freeview 2vols ./b0.nii.gz ./5tt2vis.nii.gz
python -m $SNAPSHOT --snapshot_name b0_5tt --ras_transform 2vols ./b0.nii.gz ./5tt2vis.nii.gz

mrview ./t1-in-d.nii.gz -overlay.load ./gmwmi.mif -overlay.opacity 0.3 -mode 2 &
mrconvert ./gmwmi.mif ./gmwmi.nii.gz -force
#source snapshot.sh use_freeview 2vols ./t1-in-d.nii.gz ./gmwmi.nii.gz
python -m $SNAPSHOT --snapshot_name t1_gmwmi_in_d --ras_transform 2vols ./t1-in-d.nii.gz ./gmwmi.nii.gz

#TODO: Visual checks (screenshots with freeview):


if [ "$DWI_MULTI_SHELL" = "no" ]
then
    #Estimate response function (single-tissue, single-shell)
    #???DOESN'T RECOGNIZE THIS: -nthreads 2
    dwi2response tournier ./dwi.mif ./response.txt -mask ./mask.mif -force

    #Perform spherical deconvolution to get fiber orientation distributions
    #dwi2fod ./dwi.mif ./response.txt ./wm_fod.mif -mask mask.mif -nthreads $MRTRIX_THRDS
    #dwi2fod: [ERROR] expected exactly 3 arguments (4 supplied) if I give csd in the command as below:
    dwi2fod csd ./dwi.mif ./response.txt ./wm_fod.mif -mask ./mask.mif -nthreads $MRTRIX_THRDS -force

else

    #Estimate the response functions (multi-tissue, multi-shell)
    dwi2response msmt_5tt ./dwi.mif ./5tt.mif ./RF_WM.txt ./RF_GM.txt ./RF_CSF.txt –voxels ./RF_voxels.mif -nthreads $MRTRIX_THRDS -force

    #Perform Multi-Shell, Multi-Tissue Constrained Spherical Deconvolution
    msdwi2fod msmt_csd ./dwi.mif ./RF_WM.txt ./wm_fod.mif ./RF_GM.txt ./GM.mif ./RF_CSF.txt ./CSF.mif -mask ./mask.mif -nthreads $MRTRIX_THRDS -force

    #Visual checks (interactive):
    # check appropriateness of response function #voxel selections
    mrconvert ./wm_fod.mif - -coord 3 0 | mrcat ./CSF.mif ./GM.mif - tissueRGB.mif -axis 3 -force
    #visual check:
    mrview ./tissueRGB.mif -odf.load_sh ./wm_fod.mif &

fi

#Visual check (interactive):
mrview ./t1-in-d.nii.gz -odf.load_sh ./wm_fod.mif
#TODO: Visual checks (screenshots with freeview):

popd
