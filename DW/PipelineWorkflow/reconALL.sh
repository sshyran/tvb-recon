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


#Post recon-all:

#Generate nifti files with good orientation
for vol in aparc+aseg aseg norm orig wm
do
    #Not really necessary if we have already created .nii.gz files at each previous step:
    mri_convert $MRI/$vol.mgz $MRI/$vol.nii.gz --out_orientation RAS -rt nearest
    fslreorient2std $MRI/$vol.nii.gz $MRI/$vol-reo.nii.gz
    mv $MRI/$vol-reo.nii.gz $MRI/$vol.nii.gz
done

#Get the wm surface for the cortical regions:
for h in lh rh
do
    #Transform from surface tkRAS coordinates to scanner RAS ones:
    mris_convert --to-scanner $SURF/$h.white $SURF/$h.white-ras
    #Create a white surface annotation
    cp $LABEL/$h.$DEFAULT_APARC.annot $LABEL/$h.white.annot
done

TMP=$ASEG_SURFS/tmp
if [ ! -d $TMP ]
then
    mkdir $TMP
fi
#Create surfaces around the aseg structures
for l in $ASEG_LIST
do
    lab0=$(printf %06d $l)

    # Pre-tessellate
    echo "==> Pre-tessellating: ${s}, ${lab0}"
    mri_pretess $MRI/aparc+aseg.mgz ${lab0} $MRI/norm.mgz $TMP/aseg-${lab0}.mgz

    # Tessellate
    echo "==> Tessellating: ${s}, ${lab0}"
    mri_tessellate $TMP/aseg-${lab0}.mgz ${lab0} $TMP/aseg-${lab0}-notsmooth

    mris_extract_main_component $TMP/aseg-${lab0}-notsmooth $TMP/aseg-${lab0}-notsmooth_main

    # Smooth
    echo "==> Smoothing: ${s}, ${lab0}"
    mris_smooth -nw $TMP/aseg-${lab0}-notsmooth_main $ASEG_SURFS/aseg-${lab0}

    #Transform from surface tkRAS coordinates to scanner RAS ones:
    mris_convert --to-scanner $ASEG_SURFS/aseg-${lab0} $ASEG_SURFS/aseg-ras-${lab0}

done
#rm -r $TMP


#Concatenate surfaces and create annotation:
for aseg in aseg aseg-ras
do
    python -c "import reconutils; reconutils.aseg_surf_conc_annot('$ASEG_SURFS/$aseg','$SURF/lh.$aseg','$LABEL/lh.$aseg.annot','$ASEG_LIST_LH_BS',lut_path='$FREESURFER_HOME/FreeSurferColorLUT.txt')"
    python -c "import reconutils; reconutils.aseg_surf_conc_annot('$ASEG_SURFS/$aseg','$SURF/rh.$aseg','$LABEL/rh.$aseg.annot','$ASEG_LIST_RH',lut_path='$FREESURFER_HOME/FreeSurferColorLUT.txt')"
done



