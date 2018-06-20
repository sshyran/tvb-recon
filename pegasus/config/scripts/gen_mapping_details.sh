#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}

export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh

if [ $# -eq 10 ]
then
python -m tvb.recon.qc.mapping_details -p default $@

else
python -m tvb.recon.qc.mapping_details -p $@
fi