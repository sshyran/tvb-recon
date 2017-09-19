#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}

export FREESURFER_HOME
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh
export FIGS=$PWD
export SNAPSHOT_NUMBER=$1

python -m tvb.recon.qc.snapshot ${@:2}