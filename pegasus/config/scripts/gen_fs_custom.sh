#!/usr/bin/env bash

#TODO: extract python3 env
source //anaconda/bin/activate tvb_recon_python3_env

export FREESURFER_HOME
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh
export SUBJECT=$2

python -m tvb.recon.qc.gen_fs_custom $PWD $1
