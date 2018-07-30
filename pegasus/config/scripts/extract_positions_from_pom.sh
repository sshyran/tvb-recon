#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}

python <<EOF
from tvb.recon.algo.mrielec_pos import read_write_pom_files
read_write_pom_files("$1", "$2", "$3", "$4")
EOF