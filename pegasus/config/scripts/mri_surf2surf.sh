#!/usr/bin/env bash

export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh
export TRGSUBJECT=$4
export SUBJECT=$2

f=$PWD

mri_surf2surf $@

if [ $8 == "pial" ] || [ $8 == "white" ]
then
    cd ${SUBJECTS_DIR}/${TRGSUBJECT}/surf
    if [ $6 == "lh" ]
    then
        mv lh.$8-${TRGSUBJECT} $f
    else
        mv rh.$8-${TRGSUBJECT} $f
    fi
else
    cd ${SUBJECTS_DIR}/${TRGSUBJECT}/label
    if [ $6 == "lh" ]
    then
        cp lh.aparc-${TRGSUBJECT}.annot $f
    else
        cp rh.aparc-${TRGSUBJECT}.annot $f
    fi
fi