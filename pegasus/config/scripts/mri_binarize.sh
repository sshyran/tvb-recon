#!/usr/bin/env bash

export HOME
export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh
source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}


if [ $6 == "max" ] || [ $6 == "min" ] || [ $6 == "median" ] || [ $6 == "mean" ]
then

stats=$(python <<EOF
import sys
from tvb.recon.algo.service.volume import VolumeService
stats=VolumeService().compute_vxl_stats("$2", "$6")
sys.stdout.write(str(253))
EOF)

echo mri_binarize ${@:1:5} ${stats} ${@:7}
mri_binarize ${@:1:5} ${stats} ${@:7}

else

echo mri_binarize $@
mri_binarize $@
fi