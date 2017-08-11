#!/usr/bin/env bash

# TODO: extract python3 environ
source //anaconda/bin/activate tvb_recon_python3_env

export FREESURFER_HOME
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh
export FIGS=$PWD
export SNAPSHOT_NUMBER=$1

python -m tvb.recon.qc.snapshot ${@:2}