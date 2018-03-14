#!/bin/bash

pushd $DMR


#Get volume labels:
tckmap ./$STRMLNS_SIFT_NO.tck ./tdi_ends-v$VOX.mif -vox $VOX -ends_only -template ./b0.nii.gz -force #vox: size of bin
mrconvert ./tdi_ends-v$VOX.mif ./tdi_ends-v$VOX.nii.gz -force

#Visual check (interactive)
mrview ./T1-in-d.nii.gz -overlay.load ./tdi_ends-v$VOX.mif

python -m $SNAPSHOT --ras_transform --snapshot_name tdi-v$VOX-T1-in-d 2vols $DMR/T1-in-d.nii.gz ./tdi_ends-v$VOX.nii.gz

#Label:
python -c "import tvb.recon.algo.reconutils; tvb.recon.algo.reconutils.label_vol_from_tdi('./tdi_ends-v$VOX.nii.gz','./tdi_lbl-v$VOX.nii.gz')"

python -m $SNAPSHOT --ras_transform --snapshot_name tdilbl-v$VOX-T1-in-d 2vols ./T1-in-d.nii.gz ./tdi_lbl-v$VOX.nii.gz

#Generate track counts and mean track lengths
tck2connectome -assignment_end_voxels ./$STRMLNS_SIFT_NO.tck ./tdi_lbl-v$VOX.nii.gz ./vol-counts$STRMLNS_SIFT_NO-v$VOX.csv -force
tck2connectome -assignment_end_voxels ./$STRMLNS_SIFT_NO.tck ./tdi_lbl-v$VOX.nii.gz -scale_length -stat_edge mean ./vol-mean_tract_lengths$STRMLNS_SIFT_NO-v$VOX.csv -force

#Clean up connectome of nodes without any connections, and create symmetric connectivity matrices
python -c "import tvb.recon.algo.reconutils; tvb.recon.algo.reconutils.remove_zero_connectivity_nodes('./tdi_lbl-v$VOX.nii.gz','./vol-counts$STRMLNS_SIFT_NO-v$VOX.csv',tract_length_path='./vol-mean_tract_lengths$STRMLNS_SIFT_NO-v$VOX.csv')"

python -m $SNAPSHOT --ras_transform --snapshot_name tdilbl-clean-v$VOX-T1-in-d 2vols ./T1-in-d.nii.gz ./tdi_lbl-v$VOX.nii.gz

popd
