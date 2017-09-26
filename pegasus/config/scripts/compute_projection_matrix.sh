#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}

export FREESURFER_HOME
export SUBJECTS_DIR
export SUBJECT=$4

python <<EOF
from tvb.recon.algo.reconutils import compute_projection_matrix
compute_projection_matrix("$1", "$2", "$3")
EOF