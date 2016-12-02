#!/bin/bash

#T1 input pre-processing
if [ "$T1_INPUT_FRMT" = "dicom" ]
then
    mri_convert $T1 $T1/t1_raw.nii.gz --out_orientation RAS -rt nearest
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
        mrconvert $T2 $T2/t2_raw.nii.gz
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
#source $SNAPSHOT vol_white_pial $MRI/T1-in-surf.nii.gz
python -m $SNAPSHOT --snapshot_name t1_white_pial --center_surface vol_white_pial $MRI/T1-in-surf.nii.gz

#and for aparc+aseg:
for h in lh rh
do
    #freeview -f $SURF/$h.inflated:annot=aparc \
    #    -viewport 3D -ss $FIGS/aparc-$h-annot.png
    python -m $SNAPSHOT --snapshot_name aparc_annot_$h surf_annot $SURF/$h.inflated $LABEL/$h.aparc.annot

done


