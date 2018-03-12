#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}

export FREESURFER_HOME
export SUBJECTS_DIR
export SUBJECT=$3

python <<EOF
from tvb.recon.algo.reconutils import generate_surface_zip
generate_surface_zip("$1", "$2")
EOF