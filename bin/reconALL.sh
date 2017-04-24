#!/bin/bash

#T1 input pre-processing
if [ "$T1_INPUT_FRMT" = "dicom" ]
then
    mrconvert $T1/DICOM/ $T1/t1_raw.nii.gz
    mri_convert $T1/t1_raw.nii.gz $T1/t1_raw.nii.gz --out_orientation RAS -rt nearest
    T1=$T1/t1_raw.nii.gz
fi
#ENDIF

#Freesurfer T1 processing
recon-all -s ${SUBJECT} -i $T1 -all -parallel -openmp $OPENMP_THRDS

#Additional processing if T2 and/or FLAIR is available
if [ "$T2_FLAG" = "yes" ]
then
    if [ "$T2_INPUT_FRMT" = "dicom" ]
    then
        mrconvert $T2/DICOM/ $T2/t2_raw.nii.gz
        mri_convert $T2/t2_raw.nii.gz $T2/t2_raw.nii.gz --out_orientation RAS -rt nearest
        T2=$T2/t2_raw.nii.gz
    fi
    recon-all s ${SUBJECT} -T2 $T2 -T2pial -autorecon3 -parallel -openmp $OPENMP_THRDS
fi

if [ "$FLAIR_FLAG" = "yes" ]
then
    if [ "$FLAIR_INPUT_FRMT" = "dicom" ]
    then
        mrconvert $FLAIR $FLAIR/flair_raw.nii.gz
        FLAIR=$FLAIR/flair_raw.nii.gz
    fi
    recon-all s ${SUBJECT} -FLAIR $FLAIR -FLAIRpial -autorecon3 -parallel -openmp $OPENMP_THRDS
fi

#Visual checks (screeshots) for brain, pial, white,
mri_convert $MRI/T1.mgz $MRI/T1-in-surf.nii.gz
if [ ! -d SNAPSHOTS_DIRECTORY_ENVIRON_VAR ]
then
    mkdir SNAPSHOTS_DIRECTORY_ENVIRON_VAR
fi
#source $SNAPSHOT vol_white_pial $MRI/T1-in-surf.nii.gz
python -m $SNAPSHOT --snapshot_name t1_white_pial --center_surface vol_white_pial $MRI/T1-in-surf.nii.gz

#and for aparc+aseg:
for parc in aparc.a2009s aparc.DKTatlas aparc
do
    for h in lh rh
    do
        #freeview -f $SURF/$h.inflated:annot=aparc \
        #    -viewport 3D -ss $FIGS/aparc-$h-annot.png
        python -m $SNAPSHOT --snapshot_name $parc-annot-$h surf_annot $SURF/$h.inflated $LABEL/$h.$parc.annot

    done
done

