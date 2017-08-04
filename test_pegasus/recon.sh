#!/usr/bin/env bash

export HOME=/Users/pipeline
export FREESURFER_HOME="/WORK/FS_NEW/freesurfer"
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh

f=$PWD

recon-all -all -parallel -openmp $3 -s $1 -i $2

cd $FREESURFER_HOME/subjects/$1/mri
cp T1.mgz $f
cp aparc+aseg.mgz $f
cp norm.mgz $f

cd ../surf
cp lh.pial $f
cp rh.pial $f

cd ../label
cp lh.aparc.annot $f
cp rh.aparc.annot $f