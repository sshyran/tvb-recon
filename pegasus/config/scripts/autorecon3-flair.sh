#!/usr/bin/env bash

export HOME
export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh

f=$PWD

recon-all -autorecon3 -FLAIR $2 -FLAIRpial -parallel -openmp $3 -s $1

cd ${SUBJECTS_DIR}/$1/mri
cp aparc+aseg.mgz $f