#!/bin/bash

pushd $DMR

#Tractography
#some more/alternative options?: -backtrack -crop_at_gmwmi -seed_dynamic wm_fod.mif -cutoff 0.06
opt="-seed_gmwmi ./gmwmi.mif -act ./5tt.mif -unidirectional -maxlength $STRMLNS_MAX_LEN -step $STRMLNS_STEP -nthreads $MRTRIX_THRDS"
tckgen ./wm_fod.mif ./$STRMLNS_NO.tck -number $STRMLNS_NO $opt -force

#SIFT filter
opt="-act ./5tt.mif -nthreads $MRTRIX_THRDS"
tcksift ./$STRMLNS_NO.tck ./wm_fod.mif ./$STRMLNS_SIFT_NO.tck -term_number $STRMLNS_SIFT_NO $opt -force

#Visual check (track density image -tdi)
#vox: size of bin
tckmap ./$STRMLNS_SIFT_NO.tck ./tdi_ends.mif -vox 1 -force
#Interactive:
#mrview ./t1-in-d.nii.gz -overlay.load ./tdi.mif -overlay.opacity 0.5
#Snapshot
mrconvert ./tdi_ends.mif ./tdi_ends.nii.gz -force
python -m $SNAPSHOT --snapshot_name t1_tdi_in_d --ras_transform 2vols ./t1-in-d.nii.gz ./tdi.nii.gz

popd
