#!/usr/bin/env bash

export HOME
export PATH=${MRTRIX_BIN}:${MRTRIX_SCRIPTS}:${PATH}

dwi2response msmt_5tt $1 $2 $3 $4 $5 -voxels $6 -nthreads $7