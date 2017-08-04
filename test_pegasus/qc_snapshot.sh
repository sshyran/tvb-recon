#!/usr/bin/env bash

source //anaconda/bin/activate tvb_recon_python3_env

export FREESURFER_HOME="/WORK/FS_NEW/freesurfer"
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh
export SUBJECT="TVB2PEG30"
export FIGS=$PWD
export SNAPSHOT_NUMBER=$1

python -m tvb.recon.qc.snapshot ${@:2}