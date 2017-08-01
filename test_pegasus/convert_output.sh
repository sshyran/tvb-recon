#!/usr/bin/env bash

# TODO find a better way
f=$PWD
w=$f/$1
t=$f/$2
l=$f/$3


source //anaconda/bin/activate tvb_recon_python3_env

export FREESURFER_HOME="/WORK/FS_NEW/freesurfer"
source "/WORK/FS_NEW/freesurfer/SetUpFreeSurfer.sh"

cd /WORK/BNM/tvb-recon/tvb-recon
python -m tvb.recon.qc.tvb_output $f $w $t $l $f