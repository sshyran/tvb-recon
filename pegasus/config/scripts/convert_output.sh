#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}

export FREESURFER_HOME
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh

python -m tvb.recon.qc.tvb_output -p $PWD $1 $2 $3 $PWD