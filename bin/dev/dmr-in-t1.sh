#!/bin/bash

export DMR_IN_T1=$SUBJ_DIR/dmr-in-t1

pushd $DMR_IN_T1

#Create brain mask
dwi2mask ./dwi.mif ./mask.mif -nthreads $MRTRIX_THRDS
#Extract bzero…
dwiextract ./dwi.mif ./b0.nii.gz -bzero -nthreads $MRTRIX_THRDS


#Anatomically constraint spherical deconvolution
5tt2gmwmi ./5tt.mif ./gmwmi.mif -nthreads $MRTRIX_THRDS

#Estimate response function (single-tissue, single-shell)
dwi2response tournier ./dwi.mif ./response.txt -mask ./mask.mif #???DOESN'T RECOGNIZE THIS: -nthreads 2

#Perform spherical deconvolution to get fiber orientation distributions
#dwi2fod ./dwi.mif ./response.txt ./wm_fod.mif -mask mask.mif -nthreads $MRTRIX_THRDS
#dwi2fod: [ERROR] expected exactly 3 arguments (4 supplied) if I give csd in the command as below:



dwi2fod csd ./dwi.mif ./response.txt ./wm_fod.mif -mask mask.mif -nthreads $MRTRIX_THRDS

export STRMLNS_NO=10M

opt="-seed_gmwmi ./gmwmi.mif -act ./5tt.mif -unidirectional -maxlength $STRMLNS_MAX_LEN -step $STRMLNS_STEP -nthreads $MRTRIX_THRDS"
tckgen ./wm_fod.mif ./$STRMLNS_NO.tck -number $STRMLNS_NO $opt

export STRMLNS_SIFT_NO=2M

#SIFT filter
opt="-act ./5tt.mif -nthreads $MRTRIX_THRDS"
tcksift ./$STRMLNS_NO.tck ./wm_fod.mif ./$STRMLNS_SIFT_NO.tck -term_number $STRMLNS_SIFT_NO $opt

#Visual check (track density image -tdi)
tckmap ./$STRMLNS_SIFT_NO.tck ./tdi.mif -vox 5 #vox: size of bin
mrview ./t1-in-d.nii.gz -overlay.load ./tdi.mif -overlay.opacity 0.5
