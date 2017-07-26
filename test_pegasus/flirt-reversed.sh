#!/usr/bin/env bash

export FSL_DIR="/WORK/FSL/fsl"
source ${FSL_DIR}/etc/fslconf/fsl.sh

/WORK/FSL/fsl/bin/flirt -applyxfm -in $1 -ref $2 -out $3 -init $4