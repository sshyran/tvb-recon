#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}

export FREESURFER_HOME
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh
export SUBJECT=$2

python -m tvb.recon.qc.gen_fs_custom $PWD $1
