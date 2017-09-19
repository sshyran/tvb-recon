#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}

export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh
export SUBJECT=$5

python -c "import tvb.recon.algo.reconutils; tvb.recon.algo.reconutils.aseg_surf_conc_annot('$PWD','$1','$2','$3',lut_path='$4')"