#!/usr/bin/env bash

export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh
export TRGSUBJECT=$5
export SUBJECT=$3

f=$PWD

mri_surf2surf ${@:2}

if [ $9 == "pial" ] || [ $9 == "white" ]
then
    cd ${SUBJECTS_DIR}/${TRGSUBJECT}/surf
    if [ $7 == "lh" ]
    then
        mv lh.$9-${TRGSUBJECT} $f
    else
        mv rh.$9-${TRGSUBJECT} $f
    fi
else
    cd ${SUBJECTS_DIR}/${TRGSUBJECT}/label
    if [ $7 == "lh" ]
    then
        cp lh.aparc-${TRGSUBJECT}$1.annot $f
    else
        cp rh.aparc-${TRGSUBJECT}$1.annot $f
    fi
fi