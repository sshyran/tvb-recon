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
	echo "No target subject specified via TRGSUBJECT: using fsaverage5."
	export TRGSUBJECT=fsaverage5
fi

# copy target to avoid modifying it
if [ ! -d $SUBJECTS_DIR/$TRGSUBJECT ]
then
	# try in fs home subjects
	if [ ! -d $FREESURFER_HOME/subjects/$TRGSUBJECT ]
	then
		echo "Cannot find folder for target subject $TRGSUBJECT."
		exit 1
	fi
	cp -r $FREESURFER_HOME/subjects/$TRGSUBJECT \
		$SUBJECTS_DIR/$SUBJECT-$TRGSUBJECT
else
	cp -r $SUBJECTS_DIR/$TRGSUBJECT $SUBJECTS_DIR/$SUBJECT-$TRGSUBJECT
fi

# paths to source and target folders
src=$SUBJ_DIR
trg=$SUBJECTS_DIR/$SUBJECT-$TRGSUBJECT

# map pial and white surfaces
for hemi in lh rh
do
	for sval in white pial inflated
	do
		# resamp src surf to trg
		mri_surf2surf \
			--srcsubject $SUBJECT \
			--trgsubject $SUBJECT-$TRGSUBJECT \
			--hemi $hemi \
			--sval-xyz $sval \
			--tval $sval.$SUBJECT \
			--tval-xyz $src/mri/T1.mgz
		
		# copy to src dir
		cp $trg/surf/$hemi.$sval.$SUBJECT \
		   $src/surf/$hemi.$sval.$TRGSUBJECT
	done
done

# check
# Visual check (screenshot)
#freeview -v $src/mri/T1.mgz -f $src/surf/{lh,rh}.{white,pial}.$TRGSUBJECT \
#         -viewport coronal -ss $FIGS/resamp-surfs-$SUBJECT-$TRGSUBJECT.png
source snapshot.sh vol_white_pial $src/mri/T1.nii.gz -resampled_surface_name $TRGSUBJECT

# resamp parcs
for parc in $parcs
do
	for hemi in lh rh
	do
		mri_surf2surf \
			--srcsubject $SUBJECT \
			--trgsubject $SUBJECT-$TRGSUBJECT \
			--hemi $hemi \
			--sval-annot $src/label/$hemi.$parc.annot \
			--tval $src/label/$hemi.$parc.annot.$TRGSUBJECT
	done
done

#Visual check (screenshot)
# plot resamped parcs
for parc in $parcs
do
    for hemi in lh rh
    do
        #freeview -f $src/surf/$hemi.pial.$TRGSUBJECT:annot=$src/label/$hemi.$parc.annot.$TRGSUBJECT \
        #        -viewport 3D -ss $FIGS/resamp-$hemi-pial-$parc-annot-$SUBJECT-$TRGSUBJECT.png
        mris_convert $src/surf/$hemi.pial.$TRGSUBJECT $src/surf/$hemi.pial.$TRGSUBJECT.gii
	    source snapshot.sh surf_annot $src/surf/$hemi.pial.$TRGSUBJECT.gii $src/label/$hemi.$parc.annot.$TRGSUBJECT
    done
done


# clean up
rm -r $SUBJECTS_DIR/$SUBJECT-$TRGSUBJECT
