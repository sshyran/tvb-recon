#!/bin/bash

#FLIRT co-registration of gmwmi with aseg

vol=aparc+aseg

pushd $SEGMENT

#Volume workflow:

#Get the voxels that lie at the external surface (or border) of the target volume structures:

#Using my python code:
#python -c "import reconutils; reconutils.vol_to_ext_surf_vol('$MRI/$vol.nii.gz',labels='$ASEG_LIST',hemi='lh rh',out_vol_path='./$vol-surf.nii.gz',labels_surf=None,labels_inner='0')"
#flirt -applyxfm -in ./$vol-surf.nii -ref $DMR/b0.nii.gz -out ./$vol-surf-in-d.nii.gz -init $DMR/t2d.mat -interp nearestneighbour

#Using freesurfer:
mris_calc $MRI/T1.mgz mul 0
for srf in white aseg
do
    for h in lh rh
    do
        mri_surf2vol --mkmask --hemi $h --surf $srf --identity $SUBJECT --template $MRI/T1.mgz --o ./$h.$srf-surf-mask.mgz --vtxvol ./$h.$srf-surf-map.mgz

        mris_calc ./out.mgz or ./$h.$srf-surf-mask.mgz

    done
done
mv ./out.mgz ./$vol-surf-mask.mgz
mris_calc $MRI/$vol.mgz masked ./$vol-surf-mask.mgz
mv ./out.mgz ./$vol-surf.mgz

for v in ./$vol-surf-mask ./$vol-surf
    do
        mri_convert ./$v.mgz ./$v.nii.gz --out_orientation RAS -rt nearest
        fslreorient2std ./$v.nii.gz ./$v-reo.nii.gz
        mv ./$v-reo.nii.gz ./$v.nii.gz
done


if [ "$SEGMENT_METHOD" = "tdi" ] || [ "$SEGMENT_METHOD" = "tdi+gwi" ] || [ "$SEGMENT_METHOD" = "tdi*gwi" ];
then
    pushd $TDI

    if [ -e $DMR/tdi_ends-v1.nii.gz ]
    then
        tdiends=$DMR/tdi_ends-v1.nii.gz
    else
        #Get volume labels:
        tckmap $DMR/$STRMLNS_SIFT_NO.tck ./tdi_ends.mif -vox 1 -ends_only #vox: size of bin
        mrconvert ./tdi_ends.mif ./tdi_ends.nii.gz
        rm ./tdi_ends.mif
        tdiends=./tdi_ends.nii.gz
    fi

    #Get the transform of tdi_ends in T1 space
    regopt="-dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -cost mutualinfo -interp nearestneighbour"
    flirt -in $tdiends -ref $MRI/$vol.nii.gz -omat ./tdi-to-$vol.mat -out ./tdi_ends-in-$vol.nii.gz $regopt
    #Invert the transform and get aparc+aseg to tdi_ends:
    convert_xfm -omat ./$vol-to-tdi.mat -inverse ./tdi-to-$vol.mat
    flirt -applyxfm -in $MRI/$vol.nii.gz -ref $tdiends -init ./$vol-to-tdi.mat -out ./$vol-in-tdi.nii.gz $regopt

    #...and binarize it to create a tdi mask with a threshold equal to a number of tracks
    mri_binarize --i ./tdi_ends-in-$vol.nii.gz --min $TDI_THR --o ./tdi_mask.nii.gz

    popd
fi

if [ "$SEGMENT_METHOD" = "gwi" ] || [ "$SEGMENT_METHOD" = "tdi+gwi" ] || [ "$SEGMENT_METHOD" = "tdi*gwi" ];
then
    pushd $GWI

    #Create a mask of vol's white matter
    mri_binarize --i $MRI/$vol.nii.gz --all-wm --o ./wm.nii.gz
    #is this really needed?:
    mri_convert ./wm.nii.gz ./wm.nii.gz --out_orientation RAS -rt nearest
    fslreorient2std ./wm.nii.gz ./wm-reo.nii
    mv ./wm-reo.nii.gz ./wm.nii.gz

    #gmwmi:
    #Anatomically constraint spherical deconvolution
    5ttgen fsl $MRI/T1.nii.gz ./5tt-in-T1.mif -force #if a brain mask is already applied: -premasked
    5tt2gmwmi ./5tt-in-T1.mif ./gmwmi-in-T1.mif -nthreads $MRTRIX_THRDS -force
    mrconvert ./gmwmi-in-T1.mif ./gmwmi-in-T1.nii.gz -force

    #Scale gmwmi-in-T1 in 0-256 for visualization reasons
    mris_calc ./gmwmi-in-T1.nii.gz mul 256
    mri_convert ./out.mgz ./gmwmi-in-T1-256.nii.gz --out_orientation RAS -rt nearest
    rm ./out.mgz

    #Visual checks
    #(interactive):
    #freeview -v $MRI/T1.nii ../$vol.nii.gz ./gmwmi-in-T1-256.nii:opacity=0.5
    #(screenshot):
    #freeview -v $MRI/T1.nii.gz ../$vol.nii.gz ./gmwmi-in-T1-256.nii.gz:opacity=0.5 -ss $FIGS/gmwmi-in-T1-$vol.png
    #source snapshot.sh 3vols $MRI/T1.nii.gz ../$vol.nii.gz ./gmwmi-in-T1-256.nii.gz

    #Register gmwmi with aparc+aseg
    regopt="-dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -cost mutualinfo -interp nearestneighbour"
    flirt -in ./gmwmi-in-T1.nii.gz -ref $MRI/$vol.nii.gz -omat ./gwi-to-$vol.mat -out ./gwi-in-$vol.nii.gz $regopt
    #tkregister2 --mov ./gmwmi-in-T1.nii.gz --targ ../$vol.nii.gz --reg register.dat --noedit --regheader
    #Resample gmwmi in aseg space [256 256 256]
    #mri_vol2vol --mov ./gmwmi-in-T1.nii.gz --targ ../$vol.nii.gz --o ./gmwmi-in-$vol.nii.gz --reg ./register.dat
    #Renormalize gmwmi in the [0.0, 1.0] interval
    #mris_calc ./gmwmi-in-$vol-resize.nii norm
    #mri_convert ./out.mgz ./gmwmi-in-$vol.nii.gz --out_orientation RAS -rt nearest
    #rm ./out.mgz
    #...and binarize it to create a gmwmi mask
    mri_binarize --i ./gwi-in-$vol.nii.gz --min $GWI_THR --o ./gmwmi-in-$vol-bin.mgz
    mri_convert ./gmwmi-in-$vol-bin.mgz ./gmwmi-in-$vol-bin.nii.gz --out_orientation RAS -rt nearest
    fslreorient2std ./gmwmi-in-$vol-bin.nii.gz ./gmwmi-in-$vol-bin-reo.nii.gz
    mv ./gmwmi-in-$vol-bin-reo.nii.gz ./gmwmi-in-$vol-bin.nii.gz

    #Visual checks
    #(interactive):
    #freeview -v $MRI/T1.nii.gz ../$vol.nii.gz ./gmwmi-in-$vol-bin.nii.gz:opacity=0.5
    #(screenshot):
    #freeview -v $MRI/T1.nii.gz ../$vol.nii.gz ./gmwmi-in-$vol-bin.nii.gz:opacity=0.5 -ss $FIGS/gmwmi-bin-in-$vol.png
    #source snapshot.sh 3vols $MRI/T1.nii.gz ../$vol.nii.gz ./gmwmi-in-$vol-bin.nii.gz

    #Create a mask by combining gmwmi-bin and wm masks with logical OR
    #fslmaths ./gmwmi-in-$vol-bin.nii.gz -max ./wm.nii.gz ./gwi_mask.nii.gz
    mris_calc ./gmwmi-in-$vol-bin.nii.gz or ./wm.nii.gz
    mri_convert ./out.mgz ./gwi_mask.nii.gz --out_orientation RAS -rt nearest
    rm ./out.mgz
    fslreorient2std ./gwi_mask.nii.gz ./gwi_mask-reo.nii.gz
    mv ./gwi_mask-reo.nii.gz ./gwi_mask.nii.gz

    popd
fi

if [ "$SEGMENT_METHOD" = "tdi" ]
then
    cp $TDI/tdi_mask.nii.gz ./mask-$SEGMENT_METHOD.nii.gz
elif [ "$SEGMENT_METHOD" = "tdi" ]
then
    cp $GWI/gwi_mask.nii.gz ./mask-$SEGMENT_METHOD.nii.gz
elif [ "$SEGMENT_METHOD" = "tdi+gwi" ] || [ "$SEGMENT_METHOD" = "tdi*gwi" ]
then
    if [ "$SEGMENT_METHOD" = "tdi+gwi" ]
    then
        mris_calc $TDI/tdi_mask.nii.gz or $GWI/gwi_mask.nii.gz
    else
        mris_calc $TDI/tdi_mask.nii.gz and $GWI/gwi_mask.nii.gz
    fi
    mri_convert ./out.mgz ./mask-$SEGMENT_METHOD.nii.gz --out_orientation RAS -rt nearest
    rm ./out.mgz
    fslreorient2std ./mask-$SEGMENT_METHOD.nii.gz ./mask-$SEGMENT_METHOD-reo.nii.gz
    mv ./mask-$SEGMENT_METHOD-reo.nii.gz ./mask-$SEGMENT_METHOD.nii.gz
fi

#Apply the mask to the border voxels

#Using my python code:
#python -c "import reconutils; reconutils.mask_to_vol('./$vol-surf.nii.gz','./mask-$SEGMENT_METHOD.nii.gz','./$vol-mask.nii.gz',labels='$ASEG_LIST',hemi='lh rh',vol2mask_path=None,vn=$VOL_VN,th=1,labels_mask=None,labels_nomask='0')"

#Using freesurfer:
mris_calc ./$vol-surf.nii.gz masked ./mask-$SEGMENT_METHOD.nii.gz
mv ./out.mgz ./$vol-mask.mgz
mri_convert ./$vol-mask.mgz ./$vol-mask.nii.gz --out_orientation RAS -rt nearest
fslreorient2std ./$vol-mask.nii.gz ./$vol-mask-reo.nii
mv ./$vol-mask-reo.nii.gz ./$vol-mask.nii.gz


#Get final masked surfaces:

tkregister2 --mov $MRI/T1.nii.gz --targ $MRI/T1.mgz --reg $MRI/ras2tkras.dat --noedit --regheader
mri_vol2vol --mov ./$vol-mask.nii.gz --targ $MRI/$vol.mgz --o ./$vol-mask.mgz --reg $MRI/ras2tkras.dat

#Sample the subcortical surfaces with the surviving voxels
tkregister2 --mov $MRI/T1.mgz --targ ./$vol-mask.mgz --reg $MRI/identity.dat --noedit --regheader

for h in lh rh
do
    python -c "import reconutils; reconutils.sample_vol_on_surf('$SURF/$h.aseg','./$vol-mask.nii.gz','$LABEL/$h.aseg.annot','./$h.aseg-mask',surf_ref_path='$SURF/$h.aseg-ras',out_surf_ref_path='./$h.aseg-mask-ras',ctx=None,vn=$SURF_VN)"

    #Sample mask volume on surface to create a surface mask
#mris_preproc --iv ./$vol-mask.mgz $MRI/identity.dat --out ./$h.aseg-mask.nii --hemi $h --niters 3 --target $SUBJECT --tal-xyz aseg

    #Get the actual surface using the mask:
#   python -c "import reconutils; reconutils.extract_mri_vol2subsurf('$SURF/$h.aseg','$LABEL/$h.aseg.annot','./$h.aseg-mask.nii',out_surf_path='./$h.aseg-mask2',out_annot_path='./$h.aseg-mask2.annot',ctx=None,labels=None,lut_path=os.path.join(FREESURFER_HOME,'FreeSurferColorLUT.txt'))"

done

#Sample the surface with the surviving voxels
for h in lh rh
do
    python -c "import reconutils; reconutils.sample_vol_on_surf('$SURF/$h.white','./$vol-mask.nii.gz','$LABEL/$h.aparc.annot','./$h.white-mask',surf_ref_path='$SURF/$h.white-ras',out_surf_ref_path='./$h.white-mask-ras',ctx='$h',vn=$SURF_VN)"
done



#Get the transform of tdi_lbl of the connectome specific voxel size in aparc+aseg space
#pushd $DMR

regopt="-dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -cost mutualinfo -interp nearestneighbour"
flirt -in $DMR/tdi_lbl-v$VOX.nii.gz -ref $MRI/$vol.nii.gz -omat ./tdilbl-v$VOX-to-$vol.mat -out ./tdi_lbl-v$VOX-in-$vol.nii.gz $regopt
#Invert the transform and get aparc+aseg to tdi_ends:
convert_xfm -omat ./$vol-to-tdilbl-v$VOX.mat -inverse ./tdilbl-v$VOX-to-$vol.mat
flirt -applyxfm -in $MRI/$vol.nii.gz -ref $DMR/tdi_lbl-v$VOX.nii.gz -init ./$vol-to-tdilbl-v$VOX.mat -out ./$vol-in-tdilbl-v$VOX.nii.gz $regopt
#tkregister2 --mov $MRI/$vol.nii.gz --targ ./tdi_lbl-v%VOX --reg ./$vol-to-tdilbl-v$VOX.dat --fsl ./$vol-to-tdilbl-v$VOX.mat --noedit

#tkregister2 --mov $MRI/T1.nii.gz --targ ./b0.nii.gz --reg ./t2d-fsl.dat --fsl ./t2d.mat --noedit
#tkregister2 --mov $MRI/T1.nii.gz --targ ./b0.nii.gz --reg ./t2d.dat --noedit --regheader
#mri_surf2surf --reg ./t2d.dat ./t1-in-d.nii.gz --hemi lh --sval-xyz white --tval-xyz ./t1-in-d.nii.gz --tval ./lh.white-dmr --s $SUBJECT
#mris_convert --to-scanner ./lh.white-dmr ./lh.white-dmr-ras
#popd

convert_xfm -omat ./t2dmgz.mat -inverse ./t2dmgz.mat

#python -c "import reconutils; reconutils.ijk2ijk_to_xyz2xyz('./tdi_lbl_v$vox.nii.gz','../$vol.nii.gz','./tdi-to-$vol.mat',out_path='./tdi-to-$vol-xyz.mat')"


#HOW? Distribute the rest of the vertices of that label to one of the clusters based on connectiviy and geodesic distance
#HOW? Distribute the voxels of that label to one of the sub-parcels (i.e., segmentation) base on Euclidian distance and potentially some similar connectivity contraint.

popd



