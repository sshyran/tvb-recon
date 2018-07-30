#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}

conda info

export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh

export FSL_DIR
export FSLOUTPUTTYPE=NIFTI
export PATH=${FSL_DIR}/bin:${PATH}
source ${FSL_DIR}/etc/fslconf/fsl.sh

if [ $# -eq 6 ]
then

python <<EOF
from tvb.recon.algo.mrielec_pos import transform
transform("$1", "$2", "$3", "$4", "$5", "$6", None)
EOF

else
python <<EOF
from tvb.recon.algo.mrielec_pos import transform
transform("$1", "$2", "$3", "$4", "$5", "$6", "$7")
EOF

fi
