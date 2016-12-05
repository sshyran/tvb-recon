#!/bin/bash

pushd $DMR


#Get volume labels:
tckmap ./$STRMLNS_SIFT_NO.tck ./tdi_ends-v$VOX.mif -vox $VOX -ends_only -template ./b0.nii.gz -force #vox: size of bin
mrconvert ./tdi_ends-v$VOX.mif ./tdi_ends-v$VOX.nii.gz -force

#Visual check (interactive)
mrview ./T1-in-d.nii.gz -overlay.load ./tdi_ends-v$VOX.mif
    source $CODE/snapshot.sh use_freeview 2vols $MRI/T1-in-d.nii.gz ./tdi_ends-v$VOX.nii.gz

#Label:
python -c "import reconutils; reconutils.label_vol_from_tdi('./tdi_ends-v$VOX.nii.gz','./tdi_lbl-v$VOX.nii.gz')"

#Generate track counts and mean track lengths
tck2connectome -assignment_end_voxels ./$STRMLNS_SIFT_NO.tck ./tdi_lbl-v$VOX.nii.gz ./vol-counts$STRMLNS_SIFT_NO-v$VOX.csv -force
tck2connectome -assignment_end_voxels ./$STRMLNS_SIFT_NO.tck ./tdi_lbl-v$VOX.nii.gz -scale_length -stat_edge mean ./vol-mean_tract_lengths$STRMLNS_SIFT_NO-v$VOX.csv -force

#Clean up connectome of nodes without any connections, and create symmetric connectivity matrices
python -c "import reconutils; reconutils.remove_zero_connectivity_nodes('./tdi_lbl-v$VOX.nii.gz','./vol-counts$STRMLNS_SIFT_NO-v$VOX.csv',tract_length_path='./vol-mean_tract_lengths$STRMLNS_SIFT_NO-v$VOX.csv')"


popd
