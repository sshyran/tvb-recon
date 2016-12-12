#!/bin/bash

# This script resamples anatomy annotations through another subject. The corresponding surfaces must have already been resampled and the environmental variable $TRGSUBJECT must have been modified and accordingly set by the logic in the resamp_surf.sh

# Parcellations to resample should be provided as arguments.
SUBPARCS=aparc100-geod aparc100-con+adj

# Make sure we process default ones
parcs=aparc
for default_parc in aparc.a2009s aparc.DKTatlas
do
	if [[ $parcs != *"$default_parc"* ]]
	then
		parcs="$default_parc $parcs"
	fi
done

parcs="$parcs $SUBPARCS"

# resamp parcs
for parc in $parcs
do
	for h in lh rh
	do
		mri_surf2surf \
			--srcsubject $SUBJECT \
			--trgsubject $TRGSUBJECT \
			--hemi $h \
			--sval-annot $LABEL/$h.$parc.annot \
			--tval $LABEL/$h.$parc-$TRGSUBJECT.annot
	done
done

#Visual check (screenshot)
# plot resamped parcs
for parc in $parcs
do
    for h in lh rh
    do
        #freeview -f $SURF/$h.inflated.$TRGSUBJECT:annot=$LABEL/$h.$parc-$TRGSUBJECT.annot \
        #        -viewport 3D -ss $FIGS/resamp-$h-inflated-$parc-$TRGSUBJECT-annot.png
	    python -m $SNAPSHOT --snapshot_name resamp-$parc-$TRGSUBJECT-$h surf_annot $SURF/$h.inflated-$TRGSUBJECT $LABEL/$h.$parc-$TRGSUBJECT.annot
    done
done


#TODO!: downsample aseg surf annotations!



# clean up but NOT before all annotations have been also downsampled!
#rm -r $trg
