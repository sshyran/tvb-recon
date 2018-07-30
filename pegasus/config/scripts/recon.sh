#!/usr/bin/env bash

export HOME
export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh

f=$PWD

# We assume that every case of rerunning should just resume recon-all without overwriting
# TODO: a proper management of recon-all
#if [ -d "${SUBJECTS_DIR}/$1" ]; then
#    if [ -f "${SUBJECTS_DIR}/$1/scripts/IsRunning.lh+rh" ]; then
#        rm ${SUBJECTS_DIR}/$1/scripts/IsRunning.lh+rh
#    fi
#    echo Resuming recon-all!
#    recon-all -all -no-isrunning -parallel -openmp $3 -s $1
#else
recon-all -all -parallel -openmp $3 -s $1 -i $2
#fi

cd ${SUBJECTS_DIR}/$1/mri
cp T1.mgz $f
cp aparc$4+aseg.mgz $f
cp norm.mgz $f
cp brain.mgz $f

cd ../surf
cp lh.pial $f
cp rh.pial $f
cp lh.white $f
cp rh.white $f

cd ../label
cp lh.aparc$4.annot $f
cp rh.aparc$4.annot $f