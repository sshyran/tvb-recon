#!/usr/bin/env bash

source //anaconda/bin/activate tvb_recon_python3_env

export FREESURFER_HOME="/WORK/FS_NEW/freesurfer"
source "/WORK/FS_NEW/freesurfer/SetUpFreeSurfer.sh"
export SUBJECT="TVB2PEG22"

python -c "import tvb.recon.algo.reconutils; tvb.recon.algo.reconutils.aseg_surf_conc_annot('$f','$1','$2','$3',lut_path='$4')"