#!/usr/bin/env bash

export HOME=/Users/pipeline
export FREESURFER_HOME="/WORK/FS_NEW/freesurfer"
source "/WORK/FS_NEW/freesurfer/SetUpFreeSurfer.sh"

f=$PWD

/WORK/FS_NEW/freesurfer/bin/recon-all -all -parallel -openmp $3 -s $1 -i $2

cd $FREESURFER_HOME/subjects/$1/mri
cp T1.mgz $f
cp aparc+aseg.mgz $f