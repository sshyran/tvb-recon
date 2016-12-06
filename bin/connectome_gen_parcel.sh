#!/bin/bash

pushd $DMR

# Volumes for connectome generation should be provided as arguments.
vols=aparc+aseg

for vol in $vols
do
    if [ "$vol" = "aparc+aseg" ]
    then
        #Generate labels for the default parcellation
        echo "compute label"
        labelconvert ./$vol-in-d.nii.gz $FREESURFER_HOME/FreeSurferColorLUT.txt $FS_DEFAULT ./$vol_lbl.nii.gz -force
        #older command:
        #labelconfig ./$vol-in-d.nii.gz $FS_DEFAULT ./$vol_lbl.nii.gz -lut_freesurfer $FREESURFER_HOME/FreeSurferColorLUT.txt
    else
        python -c "import reconutils; reconutils.simple_label_config('./$vol-in-d.nii.gz','./$vol_lbl.nii.gz')"
    fi

    #Generate track counts and mean track lengths for all parcellations
    assignment="-assignment_radial_search 2" #make a ball of 2 mm and look for the nearest node on the gmwgmi surface

    for metric in count meanlength
    do
        tck2connectome $STRMLNS_SIFT_NO.tck ./$vol_lbl.nii.gz $assignment ./$vol-counts$STRMLNS_SIFT_NO.csv -force
        tck2connectome $STRMLNS_SIFT_NO.tck ./$vol_lbl.nii.gz $assignment -scale_length -stat_edge mean ./$vol-mean_tract_lengths$STRMLNS_SIFT_NO.csv -force
    done
done

#TODO: quality control of the connectome.
#T1 or aparc+aseg as a background, and lines between the nodes for the network. Thickness for strength.

popd
