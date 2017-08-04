#!/usr/bin/env bash

export HOME
export FREESURFER_HOME
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh
export FSL_DIR
source ${FSL_DIR}/etc/fslconf/fsl.sh

bbregister --s $1 --mov $2 --o $3 --dti --reg $4 --lta $5 --fslmat $6
