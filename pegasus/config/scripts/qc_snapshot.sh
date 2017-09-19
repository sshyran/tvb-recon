#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}

export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh
export FIGS=$PWD
export SNAPSHOT_NUMBER=$1

python -m tvb.recon.qc.snapshot ${@:2}