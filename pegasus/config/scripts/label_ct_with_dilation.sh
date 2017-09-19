#!/usr/bin/env bash

source //anaconda/bin/activate tvb_recon_python3_env

export FREESURFER_HOME
export SUBJECTS_DIR
export SUBJECT=$4

mask="'$1'"
dilation="'$2'"
output="'$3'"

python <<EOF
from tvb.recon.algo.reconutils import label_with_dilation
label_with_dilation($mask, $dilation, $output)
EOF