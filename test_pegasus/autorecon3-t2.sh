#!/usr/bin/env bash

export HOME
export FREESURFER_HOME
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh

f=$PWD

recon-all -autorecon3 -T2 $2 -T2pial -parallel -openmp $3 -s $1

cd ${FREESURFER_HOME}/subjects/$1/mri
cp aparc+aseg.mgz $f