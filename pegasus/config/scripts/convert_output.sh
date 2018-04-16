#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}

export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh

if [ $3 == "" ]
then
    python -m tvb.recon.qc.tvb_output -p $PWD $PWD $1 $2 $PWD
else
    python -m tvb.recon.qc.tvb_output -p $PWD $PWD $1 $2 $PWD --atlas_suffix=$3
fi