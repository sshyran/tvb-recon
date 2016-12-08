#!/bin/bash

#Sub-parcellation and sub-segmentation
pushd $SUBJ_DIR

for area in $SUBAPARC_AREA
do
    for h in lh rh
    do
        python -c "import bnm.recon.algo.reconutils; bnm.recon.algo.reconutils.subparc_files('$h','aparc','aparc.sub$area',$area)"
        #or equivalently if PYTHONPATH is NOT set:
        #python bnm.recon.algo.reconutils.py subparc $h aparc aparc$area $area
    done


    #TODO: sub-segmentation of aseg
    #aseg2surf etc...

    #Create aparc+aseg
    #mri_aparc2aseg --s ${SUBJECT} --aseg aseg$area --annot aparc$area
    #for as long as we don't have aseg%area:
    mri_aparc2aseg --s ${SUBJECT} --aseg aseg --annot aparc.sub$area

    #This could be made already here:
    #Generate nifti file with good orientation
    #mri_convert $MRI/aparc$area+aseg$area.mgz $MRI/aparc$area+aseg$area.nii.gz --out_orientation RAS -rt nearest
    mri_convert $MRI/aparc.sub$area+aseg.mgz $MRI/aparc.sub$area+aseg.nii.gz --out_orientation RAS -rt nearest

    for h in lh rh
    do
        #freeview -f $SUBJ_DIR/surf/$h.inflated:annot=aparc.sub$area \
        #         -viewport 3D -ss $FIGS/subparc-$h-$area-annot-$SUBJECT.png
        mris_convert $SUBJ_DIR/surf/$h.inflated $SUBJ_DIR/surf/$h.inflated.gii
	    source snapshot.sh surf_annot $SUBJ_DIR/surf/$h.inflated.gii $SUBJ_DIR/label/$h.aparc.sub$area.annot
    done

done

popd
