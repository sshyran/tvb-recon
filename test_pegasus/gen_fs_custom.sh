#!/usr/bin/env bash

source //anaconda/bin/activate tvb_recon_python3_env

export FREESURFER_HOME="/WORK/FS_NEW/freesurfer"
source "/WORK/FS_NEW/freesurfer/SetUpFreeSurfer.sh"
export SUBJECT="TVB2PEG30"

python -m tvb.recon.qc.gen_fs_custom $PWD $1
