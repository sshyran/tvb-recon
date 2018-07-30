#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}

export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh

export FSL_DIR
export FSLOUTPUTTYPE=NIFTI
export PATH=${FSL_DIR}/bin:${PATH}
source ${FSL_DIR}/etc/fslconf/fsl.sh

python <<EOF
from tvb.recon.algo.mrielec_pos import coregister_elec_pom_and_mri
coregister_elec_pom_and_mri("$1", "$2", "$3", "$PWD", int("$4"), int("$5"))
EOF