#!/usr/bin/env bash

export HOME
export FSLDIR
export PATH=${FSLDIR}/bin:${PATH}
source ${FSLDIR}/etc/fslconf/fsl.sh
export PATH=${MRTRIX_BIN}:${MRTRIX_SCRIPTS}:${PATH}

dwipreproc $@