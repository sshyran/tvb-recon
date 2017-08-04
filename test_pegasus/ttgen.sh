#!/usr/bin/env bash

export HOME=/Users/pipeline
export PATH=/WORK/MRtrix/mrtrix3/release/bin:/WORK/MRtrix/mrtrix3/scripts:${PATH}
export FSLDIR=/WORK/FSL/fsl
export PATH=${FSLDIR}/bin:${PATH}
source ${FSLDIR}/etc/fslconf/fsl.sh

5ttgen fsl $1 $2