#!/bin/bash

pushd $DMR

# Volumes for connectome generation should be provided as arguments.
#vols=$*

#Generate labels for the default parcellation
echo "compute label"
labelconvert ./aparc+aseg-in-d.nii.gz $FREESURFER_HOME/FreeSurferColorLUT.txt $FS_DEFAULT ./aparc_lbl.nii.gz
#older command:
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

popd
