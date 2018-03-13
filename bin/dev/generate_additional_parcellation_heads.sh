#!/bin/bash

# Additional parcellations
aparc=$aparc.a2009s

# This script requires access to
# volumes: $aparc+aseg.mgz, T1.nii.gz in RAS, b0.nii.gz
# surfaces: lh.pial, rh.pial, lh.aseg, rh.aseg
# annotations: lh.$aparc.annot, rh.$aparc.annot, lh.aseg.annot, rh.aseg.annot
# color LUT file: FREESUFER_COLOR_LUT file

FREESURFER_HOME=??
COREG_USE=flirt
STRMLNS_SIFT_NO=5M #2M
SNAPSHOT=tvb.recon.qc.snapshot

DATA_FOLDER=??
SUBJECT=TVB3
SUBJECT_FOLDER=$DATA_FOLDER/$SUBJECT
OUTPUT=$SUBJECT_FOLDER/outputs
FIGS=$SUBJECT_FOLDER/figs
MRI=$OUTPUT # $SUBJECT_FOLDER/dmr
DMR=$OUTPUT # $SUBJECT_FOLDER/dmr
LABEL=$OUTPUT # $SUBJECT_FOLDER/label
SURF=$OUTPUT # $SUBJECT_FOLDER/surf

vols=$aparc+aseg

pushd $DMR

for vol in $vols
do

    #FLIRT co-registration of volumes with DWI

    #Not really necessary if we have already created .nii.gz files  with good orientation at the end recon-all:
    mri_convert $MRI/$vol.mgz $MRI/$vol.nii.gz --out_orientation RAS -rt nearest
    fslreorient2std $MRI/$vol.nii.gz $MRI/$vol-reo.nii.gz
    mv $MRI/$vol-reo.nii.gz $MRI/$vol.nii.gz

    if [ "$COREG_USE" = "flirt" ]
    then
        #Apply the transform from T1 to DWI for the volumes
        flirt -applyxfm -in $MRI/$vol.nii.gz -ref ./b0.nii.gz -out ./$vol-in-d.nii.gz -init ./t2d.mat -interp nearestneighbour
    else
        #Apply the transform from T1 to DWI for the volumes
        mri_vol2vol --mov $MRI/$vol.mgz --targ ./b0.nii.gz --o ./$vol-in-d.nii.gz --reg ./t2d.reg --nearest
    fi

    #Visual check (interactive):
    #-mode 2 is the view, you need & to allow more windows to open
    mrview ./b0.nii.gz -overlay.load ./$vol-in-d.nii.gz -overlay.opacity 0.3 -mode 2 &

    #Visual check (screenshot):
    python -m $SNAPSHOT --snapshot_name b0-$vol-in_d --ras_transform 2vols ./b0.nii.gz ./$vol-in-d.nii.gz
    #freeview -v ./T1-in-d.nii.gz ./b0.nii.gz:colormap=heat ./$vol-in-d.nii.gz:colormap=jet -ss $FIGS/t1-$vol-in-d.png
    python -m $SNAPSHOT --snapshot_name b0_t1-$vol-in_d --ras_transform 3vols ./t1-in-d.nii.gz ./b0.nii.gz ./$vol-in-d.nii.gz


    # Connectome generation

    # if [ "$vol" = "aparc+aseg" ]
    # then
    #     #Generate labels for the default parcellation
    #     echo "compute label"
    #     labelconvert ./$vol-in-d.nii.gz $FREESURFER_HOME/FreeSurferColorLUT.txt $FS_DEFAULT ./$vol_lbl.nii.gz -force
    #     #older command:
    #     #labelconfig ./$vol-in-d.nii.gz $FS_DEFAULT ./$vol_lbl.nii.gz -lut_freesurfer $FREESURFER_HOME/FreeSurferColorLUT.txt
    # else
    python -c "import tvb.recon.algo.reconutils; tvb.recon.algo.reconutils.simple_label_config('./$vol-in-d.nii.gz','./$vol_lbl.nii.gz')"
    # fi

    #Generate track counts and mean track lengths for all parcellations
    assignment="-assignment_radial_search 2" #make a ball of 2 mm and look for the nearest node on the gmwgmi surface

    for metric in count meanlength
    do
        tck2connectome $STRMLNS_SIFT_NO.tck ./$vol_lbl.nii.gz $assignment ./$vol-counts$STRMLNS_SIFT_NO.csv -force
        tck2connectome $STRMLNS_SIFT_NO.tck ./$vol_lbl.nii.gz $assignment -scale_length -stat_edge mean ./$vol-mean_tract_lengths$STRMLNS_SIFT_NO.csv -force
    done

    python -c "from tvb_recon.qc.mapping_details import compute_region_details;
               compute_region_details(atlas: AtlasSuffix,
                                      fs_color_lut: os.PathLike,
                                      t1: os.PathLike,
                                      lh_cort: os.PathLike,
                                      rh_cort: os.PathLike,
                                      lh_cort_annot: os.PathLike,
                                      rh_cort_annot: os.PathLike,
                                      lh_subcort: os.PathLike,
                                      rh_subcort: os.PathLike,
                                      lh_subcort_annot: os.PathLike,
                                      rh_subcort_annot: os.PathLike)

done

popd
