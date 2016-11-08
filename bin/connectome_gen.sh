#!/bin/bash

pushd $DMR

#IF surface parcellation generation:
if [ "$CONNECTOME_MODE" = "surf" ]
then

    # Volumes for connectome generation should be provided as arguments.
    #vols=$*

    #Generate labels for the default parcellation
    echo "compute label"
    labelconvert ./aparc+aseg-in-d.nii.gz $FREESURFER_HOME/FreeSurferColorLUT.txt $FS_DEFAULT ./aparc_lbl.nii.gz
    #labelconfig ./aparc+aseg-in-d.nii.gz $FS_DEFAULT ./aparc_lbl.nii.gz -lut_freesurfer $FREESURFER_HOME/FreeSurferColorLUT.txt

    #Generate labels for sub-parcellations (LOOP)
    for vol in $VOLS
    do
        python -c "import reconutils; reconutils.simple_label_config('./$vol-in-d.nii.gz','./$vol_lbl.nii.gz')"
    done

    #Generate track counts and mean track lengths for all parcellations
    assignment="-assignment_radial_search 2" #make a ball of 2 mm and look for the nearest node on the gmwgmi surface
    for vol in $VOLS
    do
        for metric in count meanlength
        do
        tck2connectome $STRMLNS_SIFT_NO.tck ./$vol_lbl.nii.gz $assignment ./$vol-counts$STRMLNS_SIFT_NO.csv
        tck2connectome $STRMLNS_SIFT_NO.tck ./$vol_lbl.nii.gz $assignment -scale_length -stat_edge mean ./$vol-mean_tract_lengths$STRMLNS_SIFT_NO.csv
        done
    done

else

    #Get volume labels:
    tckmap ./$STRMLNS_SIFT_NO.tck ./tdi_ends-v$VOX.mif -vox $VOX -ends_only -force #vox: size of bin
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

fi

popd
