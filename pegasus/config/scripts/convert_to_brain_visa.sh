#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}

export FREESURFER_HOME
export SUBJECTS_DIR
export SUBJECT=$4

python <<EOF
from tvb.recon.algo.reconutils import convert_fs_to_brain_visa
convert_fs_to_brain_visa("$1", "$2")
EOF