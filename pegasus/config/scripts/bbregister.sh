#!/usr/bin/env bash

export HOME
export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh
export FSL_DIR
source ${FSL_DIR}/etc/fslconf/fsl.sh

bbregister --s $1 --mov $2 --o $3 --dti --reg $4 --lta $5 --fslmat $6
