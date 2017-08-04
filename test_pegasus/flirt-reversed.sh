#!/usr/bin/env bash

export FSL_DIR="/WORK/FSL/fsl"
source ${FSL_DIR}/etc/fslconf/fsl.sh

flirt -applyxfm -in $1 -ref $2 -out $3 -init $4