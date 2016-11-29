#!/bin/bash

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
    #!!Probably not necesary anymore
    #mris_convert --to-scanner $ASEG_SURFS/aseg-${lab0} $ASEG_SURFS/aseg-ras-${lab0}

done
#rm -r $TMP


#Concatenate surfaces and create annotation:
#!!Probably not necessary for aseg-ras
for aseg in aseg
do
    python -c "import reconutils; reconutils.aseg_surf_conc_annot('$ASEG_SURFS/$aseg','$SURF/lh.$aseg','$LABEL/lh.$aseg.annot','$ASEG_LIST_LH_BS',lut_path='$FREESURFER_HOME/FreeSurferColorLUT.txt')"
    python -c "import reconutils; reconutils.aseg_surf_conc_annot('$ASEG_SURFS/$aseg','$SURF/rh.$aseg','$LABEL/rh.$aseg.annot','$ASEG_LIST_RH',lut_path='$FREESURFER_HOME/FreeSurferColorLUT.txt')"
done



