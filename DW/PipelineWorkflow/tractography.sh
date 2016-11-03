#!/bin/bash

pushd $DMR

#Tractography
#some more/alternative options?: -backtrack -crop_at_gmwmi -seed_dynamic wm_fod.mif -cutoff 0.06
opt="-seed_gmwmi ./gmwmi.mif -act ./5tt.mif -unidirectional -maxlength $STRMLNS_MAX_LEN -step $STRMLNS_STEP -nthreads $MRTRIX_THRDS"
tckgen ./wm_fod.mif ./$STRMLNS_NO.tck -number $STRMLNS_NO $opt

#SIFT filter
opt="-act ./5tt.mif -nthreads $MRTRIX_THRDS"
tcksift ./$STRMLNS_NO.tck ./wm_fod.mif ./$STRMLNS_SIFT_NO.tck -term_number $STRMLNS_SIFT_NO $opt

#Visual check (track density image -tdi)
tckmap ./$STRMLNS_SIFT_NO.tck ./tdi_ends.mif -vox 1 #vox: size of bin
#mrview ./t1-in-d.nii.gz -overlay.load ./tdi.mif -overlay.opacity 0.5
<<<<<<< HEAD
mrconvert ./tdi_ends.mif ./tdi_ends.nii.gz
source snapshot.sh 2vols ./t1-in-d.nii.gz ./tdi.nii.gz
=======
mrconvert ./tdi.mif ./tdi.nii.gz	
source snapshot.sh use_freeview 2vols ./t1-in-d.nii.gz ./tdi.nii.gz
>>>>>>> origin/DenisWorkflow

popd
