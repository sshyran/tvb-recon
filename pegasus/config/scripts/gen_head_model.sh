#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}
export FREESURFER_HOME
export SUBJECTS_DIR
export SUBJECT=$1

python <<EOF
from tvb.recon.algo.reconutils import gen_head_model
gen_head_model(decimated=$2)
EOF