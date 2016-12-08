#!/bin/bash

#FLIRT co-registration of gmwmi with aseg

pushd $ASEG_SURF

#aseg:
#Convert aseg to NIFTI with good orientation
mri_convert $MRI/aseg.mgz ./aseg.nii --out_orientation RAS -rt nearest
fslreorient2std ./aseg.nii ./aseg-reo.nii
mv ./aseg-reo.nii ./aseg.nii

#Pretess the target aseg structures
cp ./aseg.nii ./aseg_filled.nii
for l in $ASEG_LIST
do
    lbl=$(printf %03d $l)

    # Pre-tessellate
    echo "==> Pre-tessellating:"
    mri_pretess ./aseg_filled.nii ${lbl} $MRI/norm.mgz ./aseg_filled.nii
done

#Create surfaces around these border aseg volumes
export SUBASEG_VOL=$ASEG_SURF/aseg_filled
for l in $ASEG_LIST
do
    lbl=$(printf %03d $l)

    #Construct surface
    sh $CODE/aseg2srf_DP.sh -s $SUBJECT -l $l -d

    #Transform from surface tkRAS coordinates to scanner RAS ones:
    mris_convert --to-scanner $SUBASEG_VOL-${lbl} $SUBASEG_VOL-ras-${lbl}

done

#Concatenate surfaces and create annotation:
python -c "import bnm.recon.algo.reconutils; bnm.recon.algo.reconutils.aseg_surf_conc_annot('./aseg_filled','$ASEG_LIST','$FREESURFER_HOME/FreeSurferColorLUT.txt')"
python -c "import bnm.recon.algo.reconutils; bnm.recon.algo.reconutils.aseg_surf_conc_annot('./aseg_filled-ras','$ASEG_LIST','$FREESURFER_HOME/FreeSurferColorLUT.txt')"

#Get the voxels that lie at the external surface of the aseg structures:
python -c "import bnm.recon.algo.reconutils; bnm.recon.algo.reconutils.vol_to_ext_surf_vol('./aseg_filled.nii','$ASEG_LIST',out_vol_path='./aseg_filled_surf.nii',labels_surf='$ASEG_LIST',labels_inner='0')"


if [ "$ASEG_SURF_METHOD" = "gwi" ]
then
    pushd $GWI

    #Create a mask of aseg's white matter
    mri_binarize --i ../aseg_filled.nii --all-wm --o ./wm.nii
    mri_convert ./wm.nii ./wm.nii --out_orientation RAS -rt nearest
    fslreorient2std ./wm.nii ./wm-reo.nii
    mv ./wm-reo.nii ./wm.nii

    #gmwmi:
    #Anatomically constraint spherical deconvolution
    #5ttgen fsl $MRI/T1.nii.gz ./5tt-in-T1.mif -force #if a brain mask is already applied: -premasked
    5tt2gmwmi ./5tt-in-T1.mif ./gmwmi-in-T1.mif -nthreads $MRTRIX_THRDS -force
    mrconvert ./gmwmi-in-T1.mif ./gmwmi-in-T1.nii -force

    #Scale gmwmi-in-T1 in 0-256 for visualization reasons
    mris_calc ./gmwmi-in-T1.nii mul 256
    mri_convert ./out.mgz ./gmwmi-in-T1-256.nii --out_orientation RAS -rt nearest
    rm ./out.mgz

    #Visual checks
    #(interactive):
    #freeview -v $MRI/T1.nii ../aseg.nii ./gmwmi-in-T1-256.nii:opacity=0.5
    #(screenshot):
    #freeview -v $MRI/T1.nii ../aseg.nii ./gmwmi-in-T1-256.nii:opacity=0.5 -ss $FIGS/gmwmi-in-T1-aseg.png

    #Register gmwmi with aseg
    tkregister2 --mov ./gmwmi-in-T1.nii --targ ../aseg.nii --reg register.dat --noedit --regheader
    #Resample gmwmi in aseg space [256 256 256]
    mri_vol2vol --mov ./gmwmi-in-T1.nii --targ ../aseg.nii --o ./gmwmi-in-T1-resize.nii --reg ./register.dat
    #Renormalize gmwmi in the [0.0, 1.0] interval
    mris_calc ./gmwmi-in-T1-resize.nii norm
    #...and binarize it to create a gmwmi mask
    mri_binarize --i ./out.mgz --min $ASEG_GWI_THR --o ./gmwmi-in-T1-bin.mgz
    rm ./out.mgz
    mri_convert ./gmwmi-in-T1-bin.mgz ./gmwmi-in-T1-bin.nii --out_orientation RAS -rt nearest
    fslreorient2std ./gmwmi-in-T1-bin.nii ./gmwmi-in-T1-bin-reo.nii
    mv ./gmwmi-in-T1-bin-reo.nii ./gmwmi-in-T1-bin.nii

    #Visual checks
    #(interactive):
    #freeview -v $MRI/T1.nii ../aseg.nii ./gmwmi-in-T1-bin.nii.nii:opacity=0.5
    #(screenshot):
    #freeview -v $MRI/T1.nii ../aseg.nii ./gmwmi-in-T1-bin.nii:opacity=0.5 -ss $FIGS/gmwmi-bin-in-T1-aseg.png

    #Create a mask by combining gmwmi-bin and wm masks with logical OR
    #fslmaths ./gmwmi-in-T1-bin.nii -max ./wm.nii ./gwi_mask.nii.gz
    mris_calc ./gmwmi-in-T1-bin.nii or ./wm.nii
    mri_convert ./out.mgz ./gwi_mask.nii --out_orientation RAS -rt nearest
    rm ./out.mgz
    fslreorient2std ./gwi_mask.nii ./gwi_mask-reo.nii
    mv ./gwi_mask-reo.nii ./gwi_mask.nii

    #Get the interface of the aseg surface voxels with white matter
    python -c "import bnm.recon.algo.reconutils; bnm.recon.algo.reconutils.mask_to_vol('../aseg_filled_surf.nii','./gwi_mask.nii','../aseg_filled_surf_con_gwi.nii','$ASEG_LIST',vol2maskTrnsfrmPath=None,vn=$ASEG_GWI_VN,th=1,labels_mask='$ASEG_LIST',labels_nomask='0')"

    popd

elif [ "$ASEG_SURF_METHOD" = "tdi" ]
then
    pushd $TDI

    #tdi: rejected due to bad registrationg from DWI to T1...
    #tckmap $DMR/$STRMLNS_SIFT_NO.tck ./tdi.mif -vox 1 #vox: size of bin
    #mrconvert ./tdi.mif ./tdi.nii -vox 1,1,1
    #regopt="-dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -cost mutualinfo"
    #flirt -in ./tdi.nii -ref $MRI/T1.nii.gz -omat ./tdi2t.mat -out ./td1-in-t1.nii.gz $regopt

    #Visual check (track density image -tdi)
    tckmap $DMR/$STRMLNS_SIFT_NO.tck ./tdi.mif -vox 1 #vox: size of bin
    #mrview $DMR/t1-in-d.nii.gz -overlay.load ./tdi.mif -overlay.opacity 0.5
    mrconvert ./tdi.mif ./tdi.nii

    #...and binarize it to create a tdi mask with a threshold equal to a number of tracks
    mri_binarize --i ./tdi.nii --min $ASEG_TDI_THR --o ./tdi_mask.nii

    #Get the interface of the aseg surface voxels with tracks' ends
    python -c "import bnm.recon.algo.reconutils; bnm.recon.algo.reconutils.mask_to_vol('../aseg_filled_surf.nii','./tdi_mask.nii','../aseg_filled_surf_con_tdi.nii','$ASEG_LIST',vol2maskTrnsfrmPath='$DMR/t2d.mat',vn=$ASEG_TDI_VN,th=1,labels_mask='$ASEG_LIST',labels_nomask='0')"

    popd
fi

#Sample the aseg surfs with their gmwmi/tracts interface voxels
ASEG_CON="./aseg_filled_surf_con*.nii"
#ASEG_CON=( $ASEG_CON )
#ASEG_CON=${ASEG_CON[0]}
#echo ASEG_CON=ASEG_CON
#ASEG_CON_SURF=${ASEG_CON%.nii}
for aseg in $ASEG_CON
do
    ASEG_CON_SURF=${aseg%.nii}
    echo $ASEG_CON_SURF
    python -c "import bnm.recon.algo.reconutils; bnm.recon.algo.reconutils.sample_vol_on_surf('./aseg_filled','$aseg','$ASEG_CON_SURF','$ASEG_LIST', vn=$SURF_VN)"

    #Concatenate surfaces and create annotation:
    python -c "import bnm.recon.algo.reconutils; bnm.recon.algo.reconutils.aseg_surf_conc_annot('$ASEG_CON_SURF','$ASEG_LIST','$FREESURFER_HOME/FreeSurferColorLUT.txt')"
done

#Move individual aseg surfs files to get a more tidy directory
if [ ! -d ./aseg_filled_surfs ]
then
    mkdir ./aseg_filled_surfs
fi
mv ./aseg_filled-0* ./aseg_filled_surfs
mv ./aseg_filled-ras-0* ./aseg_filled_surfs

if [ -d ./gwi ]
then
    if [ ! -d ./gwi/surfs ]
    then
        mkdir ./gwi/surfs
        mv ./aseg_filled_surf_con_gwi-* ./gwi/surfs
    fi
fi
if [ -d ./tdi ]
then
    if [ ! -d ./tdi/surfs ]
    then
        mkdir ./tdi/surfs
        mv ./aseg_filled_surf_con_tdi-* ./tdi/surfs
    fi
fi

#And/Or create new surfaces around the interface voxels with white matter
#export SUBASEG_VOL=$ASEG_SURF/aseg_filled_surf_con
#Create surfaces around these border aseg volumes
#for l in $ASEG_LIST
#do
    #lbl=$(printf %03d $l)

    #Construct surface
    #sh $CODE/aseg2srf_DP.sh -s $SUBJECT -l $l -d
#done

#Calculate aseg surface areas

#Sub-parcel the aseg surfaces

#Sub-segment aseg by assigning voxels according to distance to the sub-parcels

#Visual check (screenshot):
#freeview -f ./aseg*:color=grey ./gw_mask.nii:opacity=0.5 -viewport 3d -ss $FIGS/gw-aseg-surfs.png

popd



