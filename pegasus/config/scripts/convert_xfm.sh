#!/usr/bin/env bash

export FSL_DIR
export PATH=${FSL_DIR}/bin:${PATH}
source ${FSL_DIR}/etc/fslconf/fsl.sh

convert_xfm $@
