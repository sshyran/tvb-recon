#!/usr/bin/env bash

export MNE_ROOT
source ${MNE_ROOT}/bin/mne_setup_sh
export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh
export SUBJECT=$1

f=$PWD

mne_watershed_bem --overwrite

cd ${SUBJECTS_DIR}/${SUBJECT}/bem/watershed

cp ${SUBJECT}_brain_surface $f
cp ${SUBJECT}_inner_skull_surface $f
cp ${SUBJECT}_outer_skin_surface $f
cp ${SUBJECT}_outer_skull_surface $f