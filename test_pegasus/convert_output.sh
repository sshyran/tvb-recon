#!/usr/bin/env bash

source //anaconda/bin/activate tvb_recon_python3_env

export FREESURFER_HOME
source ${FREESURFER_HOME}/SetUpFreeSurfer.sh

python -m tvb.recon.qc.tvb_output -p $PWD $1 $2 $3 $PWD