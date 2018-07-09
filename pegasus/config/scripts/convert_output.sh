#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}

export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh

if [ $# -eq 3 ]
then
python -m tvb.recon.qc.tvb_output -p $3 $PWD $PWD $1 $2 $PWD

else
python -m tvb.recon.qc.tvb_output -p default $PWD $PWD $1 $2 $PWD

fi