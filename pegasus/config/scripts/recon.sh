#!/usr/bin/env bash

export HOME
export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh

f=$PWD

recon-all -all -parallel -openmp $3 -s $1 -i $2

cd ${SUBJECTS_DIR}/$1/mri
cp T1.mgz $f
cp aparc+aseg.mgz $f
cp norm.mgz $f
cp brain.mgz $f

cd ../surf
cp lh.pial $f
cp rh.pial $f

cd ../label
cp lh.aparc.annot $f
cp rh.aparc.annot $f