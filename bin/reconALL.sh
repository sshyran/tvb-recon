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
source $CODE/snapshot.sh vol_white_pial $MRI/T1-in-surf.nii.gz

#and for aparc+aseg:
for h in lh rh
do
    #freeview -f $SUBJ_DIR/surf/$h.inflated:annot=aparc \
    #    -viewport 3D -ss $FIGS/aparc-$h-aparc-annot-$SUBJECT.png
    mris_convert $SUBJ_DIR/surf/$h.inflated $SUBJ_DIR/surf/$h.inflated.gii
    source $CODE/snapshot.sh surf_annot $SUBJ_DIR/surf/$h.inflated.gii $SUBJ_DIR/label/$h.aparc.annot
done


