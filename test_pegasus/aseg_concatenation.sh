#!/usr/bin/env bash

# TODO find a better way
f=$PWD
s=$f/$1
a=$f/$2
l=$f/$4


source //anaconda/bin/activate tvb_recon_python3_env

export FREESURFER_HOME="/WORK/FS_NEW/freesurfer"
source "/WORK/FS_NEW/freesurfer/SetUpFreeSurfer.sh"
export SUBJECT="TVB2PEG22"

cd /WORK/BNM/tvb-recon/tvb-recon

python -c "import tvb.recon.algo.reconutils; tvb.recon.algo.reconutils.aseg_surf_conc_annot('$f','$s','$a','$3',lut_path='$l')"