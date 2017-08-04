#!/usr/bin/env bash

export HOME
export PATH=${MRTRIX_BIN}:${MRTRIX_SCRIPTS}:${PATH}
export FSLDIR
export PATH=${FSLDIR}/bin:${PATH}
source ${FSLDIR}/etc/fslconf/fsl.sh

5ttgen fsl $1 $2