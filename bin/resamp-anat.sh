#!/bin/bash

# This script resamples anatomy through another subject. By default,
# fsaverage5 subject is used, which results in a lower resolution
# cortical surface.

# Parcellations to resample should be provided as arguments.
parcs=aparc
parcs="$SUBPARCS $parcs"

# Make sure we process default ones
for default_parc in aparc.a2009s aparc.DKTatlas
do
	if [[ $parcs != *"$default_parc"* ]]
	then
		parcs="$default_parc $parcs"
	fi
done

# establish target anatomy
if [ -z "$TRGSUBJECT" ]
then
	echo "No target subject specified via TRGSUBJECT. Creating a new one with decimation factor:"
    echo $DECIM_FACTOR
    export TRGSUBJECT=$SUBJECT-RESAMP_$DECIM_FACTOR
    echo "TRGSUBJECT="
    echo $TRGSUBJECT
fi


#If this target subject doesn't exist...
if [ ! -d $SUBJECTS_DIR/$TRGSUBJECT ]
then
    #...and if it is not one of the default freesurfer subjects
    if [ ! -d $FREESURFER_HOME/subjects/$TRGSUBJECT ]
    then
        #Create necessary folders for the resampled target subject
        trg=$SUBJECTS_DIR/$TRGSUBJECT
        mkdir $trg
        trgsurf=$trg/surf
        mkdir $trgsurf
        trglbl=$trg/label
        mkdir $trglbl
        #Decimate the original l/rh.sphere.reg
        for h in lh rh
        do
            mris_decimate -d $DECIM_FACTOR $SURF/$h.sphere.reg $trgsurf/$h.sphere.reg
        done
    else
        #We can use the one of the default subject
        cp -r $FREESURFER_HOME/subjects/$TRGSUBJECT $SUBJECTS_DIR/$SUBJECT-$TRGSUBJECT
        export TRGSUBJECT=$SUBJECT-$TRGSUBJECT
        echo "TRGSUBJECT="
        echo TRGSUBJECT
        #paths to target folders
        trg=$SUBJECTS_DIR/$TRGSUBJECT
        trgsurf=$trg/surf
        trglbl=$trg/label
    fi
else
    #If the target subject exists:
    #If we are inside the freesurfer home folder, we need to make a copy to protect the original
    if [ $SUBJECTS_DIR = $FREESURFER_HOME/subjects ]
    then
        cp -r $FREESURFER_HOME/subjects/$TRGSUBJECT $SUBJECTS_DIR/$SUBJECT-$TRGSUBJECT§
        export TRGSUBJECT=$SUBJECT-$TRGSUBJECT
        echo "TRGSUBJECT="
        echo TRGSUBJECT
    fi
    #paths to target folders
    trg=$SUBJECTS_DIR/$TRGSUBJECT
    trgsurf=$trg/surf
    trglbl=$trg/label
fi

# paths to source folders
src=$SUBJ_DIR

# map pial and white surfaces
for hemi in lh rh
do
	for sval in white pial inflated
	do
		# resamp src surf to trg
		mri_surf2surf \
			--srcsubject $SUBJECT \
			--trgsubject $TRGSUBJECT \
			--hemi $hemi \
			--sval-xyz $sval \
			--tval $sval-$SUBJECT \
			--tval-xyz $MRI/T1.mgz
		
		# copy to src dir
		cp $trgsurf/$hemi.$sval-$SUBJECT \
		   $SURF/$hemi.$sval-$TRGSUBJECT
	done
done

# check
# Visual check (screenshot)
#freeview -v $MRI/T1.mgz -f $SURF/{lh,rh}.{white,pial}-$TRGSUBJECT \
#         -viewport coronal -ss $FIGS/resamp-surfs-$TRGSUBJECT.png
python -m $SNAPSHOT --center_surface --snapshot_name resamp_white_pial_t1_$TRGSUBJECT vol_white_pial $MRI/T1.nii.gz -resampled_surface_name $TRGSUBJECT

# resamp parcs
for parc in $parcs
do
	for hemi in lh rh
	do
		mri_surf2surf \
			--srcsubject $SUBJECT \
			--trgsubject $TRGSUBJECT \
			--hemi $hemi \
			--sval-annot $LABEL/$hemi.$parc.annot \
			--tval $LABEL/$hemi.$parc-$TRGSUBJECT.annot
	done
done

#Visual check (screenshot)
# plot resamped parcs
for parc in $parcs
do
    for hemi in lh rh
    do
        #freeview -f $SURF/$hemi.inflated.$TRGSUBJECT:annot=$LABEL/$hemi.$parc-$TRGSUBJECT.annot \
        #        -viewport 3D -ss $FIGS/resamp-$hemi-inflated-$parc-$TRGSUBJECT-annot.png
	    python -m $SNAPSHOT --snapshot_name resamp-$parc-$TRGSUBJECT surf_annot $SURF/$hemi.inflated-$TRGSUBJECT $LABEL/$hemi.$parc-$TRGSUBJECT.annot
    done
done

#Downsample now the aseg surfs
for h in rh lh
do
    mris_decimate -d $DECIM_FACTOR $SURF/$h.aseg $SURF/$h.aseg-$TRGSUBJECT
done
#Screenshot:
python -m $SNAPSHOT --center_surface --snapshot_name resamp_aseg_t1_$TRGSUBJECT vol_surf $MRI/T1.nii.gz $SURF/{lh,rh}.aseg $SURF/{lh,rh}.aseg-$TRGSUBJECT

#TODO!: downsample aseg surf annotations!



# clean up
rm -r $trg
