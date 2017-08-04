#!/usr/bin/env bash

export HOME=/Users/pipeline
export FREESURFER_HOME="/WORK/FS_NEW/freesurfer"
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh
export FSL_DIR="/WORK/FSL/fsl"
source ${FSL_DIR}/etc/fslconf/fsl.sh

bbregister --s $1 --mov $2 --o $3 --dti --reg $4 --lta $5 --fslmat $6
