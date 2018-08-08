#!/usr/bin/env bash

export HOME
export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh

mri_binarize $@