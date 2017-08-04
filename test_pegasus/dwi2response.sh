#!/usr/bin/env bash

export HOME
export PATH=${MRTRIX_BIN}:${MRTRIX_SCRIPTS}:${PATH}

dwi2response tournier $1 $2 -mask $3