#!/usr/bin/env bash

export HOME=/Users/pipeline
export FREESURFER_HOME="/WORK/FS_NEW/freesurfer"
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh

f=$PWD

recon-all -autorecon3 -FLAIR $2 -FLAIRpial -parallel -openmp $3 -s $1

cd ${FREESURFER_HOME}/subjects/$1/mri
cp aparc+aseg.mgz $f