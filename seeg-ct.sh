seeg_dir=$SUBJECTS_DIR/$SUBJECT/seeg
mkdir $seeg_dir
pushd $seeg_dir

mri_convert $elec_image elec.nii.gz

for image in T1 brain aparc+aseg
do
	mri_convert ../mri/$image.mgz $image.nii.gz
done

for image in T1 elec brain aparc+aseg
do
	fslreorient2std $image.nii.gz $image-reo.nii.gz
done

regopts="-dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -cost mutualinfo"
flirt -in elec-reo.nii.gz -ref T1-reo.nii.gz -omat elec_to_t1.mat -out elec_in_t1.nii.gz $regopts

convert_xfm -omat elec_to_t1_inv.mat -inverse elec_to_t1.mat 

for image in T1 brain aparc+aseg
do
	flirt -applyxfm -in $image-reo.nii.gz -ref elec-reo.nii -out $image-in-elec.nii.gz -init elec_to_t1_inv.mat 
	# if aparc+aseg, need  -interp nearestneighbour
done


mri_binarize --erode 8 --i brain-in-elec.nii.gz --o eroded-mask.nii.gz --min 10
# load freeview, using T1 as background, elec w/ erorded as mask, it fits fine

# now binarize elec using eroded mask
mri_binarize --i elec_reo.nii --o elec_bin.nii --min 100 --mask boo.nii.gz

mri_binarize --i elec-reo.nii --o elec-mask-bin.nii --min 100 --mask eroded-mask.nii.gz

# works well, though we can't reliable distinguish contacts with our CTs or MRIs, so
# we dilate
mri_binarize --i elec_reo.nii --o elec_bin_d4.nii --min 500 --dilate 4 --erode 2 --mask boo.nii.gz

# this works nicely, now we need to fit a linear geometry to the blobs

# could do a graph cut, then svd coords of each
# scipy.ndimage.label works for the cut, but electrodes get too close sometimes
# second singular val should tell whether two elecs together or not

# dilated is a little bloated, maybe can use only to label voxels of elec-mask-bin?
# is just array mult in python log

# project label voxels onto first PC, min max gives extent of electrode
# 

# TODO threshold for MRI electrodes won't be valid
