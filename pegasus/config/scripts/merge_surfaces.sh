#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}

export FREESURFER_HOME
export SUBJECTS_DIR
export SUBJECT=$4

python <<EOF
from tvb.recon.algo.reconutils import merge_surfs
merge_surfs("$1", "$2", "$3")
EOF
