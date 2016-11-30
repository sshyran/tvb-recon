#!/bin/bash

#FLIRT co-registration of gmwmi with aseg

vol=aparc+aseg

pushd $SEGMENT


#Volume workflow:


if [ "$SEGMENT_METHOD" = "tdi" ] || [ "$SEGMENT_METHOD" = "tdi+gwi" ] || [ "$SEGMENT_METHOD" = "tdi*gwi" ];
then
    pushd $TDI

    if [ -e $DMR/tdi_ends-v1.nii.gz ]
    then
        tdiends=$DMR/tdi_ends-v1.nii.gz
    else
        #Get volume labels:
        tckmap $DMR/$STRMLNS_SIFT_NO.tck ./tdi_ends.mif -vox 1.0 -ends_only -template $DMR/b0.nii.gz
        mrconvert ./tdi_ends.mif ./tdi_ends.nii.gz
        tdiends=./tdi_ends.nii.gz
    fi

    #Get the transform of tdi_ends in T1 space
    #QUESTION! Can we afford the loss of accuracy due to volume resampling from diffusion space and resolution to those of T1?

    if [ "$COREG_USE" = "flirt" ]
    then
        flirt -applyxfm -in $tdiends -ref $MRI/T1.nii.gz -init $DMR/t2d.mat -out ./tdi_ends-in-t1.nii.gz
    else
        mri_vol2vol --mov $tdiends --targ $MRI/T1.mgz --o ./tdi_ends-in-t1.mgz --reg ./d2t.reg
        mri_convert ./tdi_ends-in-t1.mgz ./tdi_ends-in-t1.nii.gz -out_orientation ras
        rm ./tdi_ends-in-t1.mgz
    fi

    #...and binarize it to create a tdi mask with a threshold equal to the number of tracks
    mri_binarize --i ./tdi_ends-in-t1.nii.gz --min $TDI_THR --o ./tdi_mask.nii.gz

    popd
fi

if [ "$SEGMENT_METHOD" = "gwi" ] || [ "$SEGMENT_METHOD" = "tdi+gwi" ] || [ "$SEGMENT_METHOD" = "tdi*gwi" ];
then
    pushd $GWI

    #Create a mask of vol's white matter
    mri_binarize --i $MRI/$vol.nii.gz --all-wm --o ./wm.nii.gz
    #is this really needed?:
    mri_convert ./wm.nii.gz ./wm.nii.gz --out_orientation RAS -rt nearest
    ##!!Probably not necesary anymore
    #fslreorient2std ./wm.nii.gz ./wm-reo.nii
    #mv ./wm-reo.nii.gz ./wm.nii.gz

    #gmwmi:
    #Anatomically constraint spherical deconvolution
    #5ttgen fsl $MRI/T1.nii.gz ./5tt-in-t1.mif -force #if a brain mask is already applied: -premasked
    5tt2gmwmi ./5tt-in-t1.mif ./gmwmi-in-t1.mif -nthreads $MRTRIX_THRDS -force
    mrconvert ./gmwmi-in-t1.mif ./gmwmi-in-t1.nii.gz -force

    #Scale gmwmi-in-T1 in 0-256 for visualization reasons
    mris_calc ./gmwmi-in-t1.nii.gz mul 256
    mri_convert ./out.mgz ./gmwmi-in-t1-256.nii.gz --out_orientation RAS -rt nearest
    rm ./out.mgz

    #Visual checks
    #(interactive):
    #freeview -v $MRI/T1.nii ../$vol.nii.gz ./gmwmi-in-t1-256.nii:opacity=0.5
    #(screenshot):
    #freeview -v $MRI/T1.nii.gz ../$vol.nii.gz ./gmwmi-in-t1-256.nii.gz:opacity=0.5 -ss $FIGS/gmwmi-in-t1-$vol.png
    #source snapshot.sh 3vols $MRI/T1.nii.gz ../$vol.nii.gz ./gmwmi-in-t1-256.nii.gz

    #Resample and register gmwmi with aparc+aseg
    #tkregister2 --mov ./gmwmi-in-t1.nii.gz --targ $MRI/T1.nii.gz --reg ./resamp_gw-in-t1.dat --noedit --regheader
    #Resample gmwmi in aseg space [256 256 256]
    #mri_vol2vol --mov ./gmwmi-in-t1.nii.gz --targ $MRI/T1.nii.gz --o ./gmwmi-in-t1-resamp.nii.gz --reg ./resamp_gw-in-t1.dat
    mri_vol2vol --mov ./gmwmi-in-t1.nii.gz --targ $MRI/T1.nii.gz --o ./gmwmi-in-t1-resamp.nii.gz --regheader
    #Renormalize gmwmi in the [0.0, 1.0] interval
    mris_calc ./gmwmi-in-t1-resamp.nii.gz norm
    mri_convert ./out.mgz ./gmwmi-in-t1-resamp-norm.nii.gz --out_orientation RAS -rt nearest
    #rm ./out.mgz
    #...and binarize it to create a gmwmi mask
    mri_binarize --i ./gmwmi-in-t1-resamp-norm.nii.gz --min $GWI_THR --o ./gmwmi-in-t1-bin.mgz
    mri_convert ./gmwmi-in-t1-bin.mgz ./gmwmi-in-t1-bin.nii.gz --out_orientation RAS -rt nearest
    ##!!Probably not necesary anymore
    #fslreorient2std ./gmwmi-in-t1-bin.nii.gz ./gmwmi-in-t1-bin-reo.nii.gz
    #mv ./gmwmi-in-t1-bin-reo.nii.gz ./gmwmi-in-t1-bin.nii.gz

    #Visual checks
    #(interactive):
    #freeview -v $MRI/T1.nii.gz ../$vol.nii.gz ./gmwmi-in-t1-bin.nii.gz:opacity=0.5
    #(screenshot):
    #freeview -v $MRI/T1.nii.gz ../$vol.nii.gz ./gmwmi-in-t1-bin.nii.gz:opacity=0.5 -ss $FIGS/gmwmi-bin-in-$vol.png
    #source snapshot.sh 3vols $MRI/T1.nii.gz ../$vol.nii.gz ./gmwmi-in-t1-bin.nii.gz

    #Create a mask by combining gmwmi-bin and wm masks with logical OR
    mris_calc ./gmwmi-in-t1-bin.nii.gz or ./wm.nii.gz
    mri_convert ./out.mgz ./gwi_mask.nii.gz --out_orientation RAS -rt nearest
    rm ./out.mgz
    ##!!Probably not necesary anymore
    #fslreorient2std ./gwi_mask.nii.gz ./gwi_mask-reo.nii.gz
#   mv ./gwi_mask-reo.nii.gz ./gwi_mask.nii.gz

    popd
fi

#Create the final mask

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
    ##!!Probably not necesary anymore
    #fslreorient2std ./mask-$SEGMENT_METHOD.nii.gz ./mask-$SEGMENT_METHOD-reo.nii.gz
    #mv ./mask-$SEGMENT_METHOD-reo.nii.gz ./mask-$SEGMENT_METHOD.nii.gz
fi


#If segmentation is performed via parcellation of the surface, maybe it is not necessary to work on the surface voxels of aparc+aseg volume. In this case we directly mask the aparc+aseg.
if [ "$APARC_SURF" = "yes" ]
then
    #Apply the mask to the border voxels

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
        ##!!Probably not necesary anymore
        #fslreorient2std ./$v.nii.gz ./$v-reo.nii.gz
        #mv ./$v-reo.nii.gz ./$v.nii.gz
    done

    mask_this_vol=$vol-surf
else
    #Apply the mask to the whole aparc+aseg, not only to the border voxels
    mask_this_vol=$vol
done

#Using my python code:
#python -c "import reconutils; reconutils.mask_to_vol('./$vol-surf.nii.gz','./mask-$SEGMENT_METHOD.nii.gz','./$vol-mask.nii.gz',labels='$ASEG_LIST',hemi='lh rh',vol2mask_path=None,vn=$VOL_VN,th=1,labels_mask=None,labels_nomask='0')"

#Using freesurfer:
mris_calc ./$mask_this_vol.nii.gz masked ./mask-$SEGMENT_METHOD.nii.gz
mri_convert ./out.mgz ./$vol-mask.nii.gz --out_orientation RAS -rt nearest
rm ./out.mgz
##!!Probably not necesary anymore
#fslreorient2std ./$vol-mask.nii.gz ./$vol-mask-reo.nii
#mv ./$vol-mask-reo.nii.gz ./$vol-mask.nii.gz


#Get final masked surfaces by sampling the surface with the surviving voxels:
#Cortical:
#Give an empty list for add_lbl, if you want cerebral white matter to be masked out
for h in lh rh
do
    python -c "import reconutils; reconutils.sample_vol_on_surf('$SURF/$h.white','./$vol-mask.nii.gz','$LABEL/$h.aparc.annot','./$h.white-mask','$vox2rastkr_path',ctx='$h',vn=$SURF_VN,add_lbl=[2,41])"
done
#Sub-cortical:
for h in lh rh
do
    python -c "import reconutils; reconutils.sample_vol_on_surf('$SURF/$h.aseg','./$vol-mask.nii.gz','$LABEL/$h.aseg.annot','./$h.aseg-mask','$vox2rastkr_path',ctx=None,vn=$SURF_VN,add_lbl=[])"
done


#If connectivity similarity is one of the criteria for parcellation
if [[ $SUBAPARC_MODE == *"con"* ]]
then

    out_consim_path=./consim-vol-counts$STRMLNS_SIFT_NO-v$VOX.npy
    #Compute the voxel connectivity similarity
    #This takes a lot of time...
    python -c "import reconutils; reconutils.node_connectivity_metric('$DMR/vol-counts$STRMLNS_SIFT_NO-v$VOX.npy',metric='cosine', mode='sim',out_consim_path='$out_consim_path')"

    #Sample the connectome volume to T1 space
    if [ "$COREG_USE" = "flirt" ]
    then
        flirt -applyxfm -in $DMR/tdi_lbl-v$VOX.nii.gz -ref $MRI/T1.nii.gz -init $DMR/t2d.mat -out ./tdi_lbl-v$VOX-in-t1.nii.gz -interp nearestneighbour
        ref_vol_path=./tdi_lbl-v$VOX-in-t1.nii.gz
    else
        mri_vol2vol --mov $DMR/tdi_lbl-v$VOX.nii.gz --targ $MRI/T1.mgz --o ./tdi_lbl-v$VOX-in-t1.mgz --reg ./d2t.reg --nearest
        mri_convert ./tdi_lbl-v$VOX-in-t1.mgz ./tdi_lbl-v$VOX-in-t1.nii.gz-out_orientation ras
        ref_vol_path=./tdi_lbl-v$VOX-in-t1.mgz

    fi

else
    ref_vol_path=''
    vox2rastkr_path=
fi

for h in lh rh
do
    python -c "import reconutils; reconutils.connectivity_geodesic_subparc('$SURF/$h.white','$LABEL/$h.aparc.annot','./$h.white-mask-idx.npy',out_annot_path='$SURF/$h.aparc$SUBAPARC_AREA-$SUBAPARC_MODE.annot', ref_vol_path='$ref_vol_path',consim_path='$out_consim_path',parc_area=$SUBAPARC_AREA, labels=None,hemi='$h', mode='$SUBAPARC_MODE', vox2rastkr_path='$vox2rastkr_path',lut_path=os.path.join('$FREESURFER_HOME','FreeSurferColorLUT.txt'))"
done

for h in lh rh
do
    python -c "import reconutils; reconutils.connectivity_geodesic_subparc('$SURF/$h.aseg','$LABEL/$h.aseg.annot','./$h.white-mask-idx.npy',out_annot_path='$SURF/$h.aseg$SUBAPARC_AREA-$SUBAPARC_MODE.annot', ref_vol_path='$ref_vol_path',consim_path='$out_consim_path',parc_area=$SUBAPARC_AREA, labels='${ASEG_LIST_$h}',hemi=None, mode='$SUBAPARC_MODE', vox2rastkr_path='$vox2rastkr_path',lut_path=os.path.join('$FREESURFER_HOME','FreeSurferColorLUT.txt'))"
done
popd



