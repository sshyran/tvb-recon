#!/usr/bin/env bash

source ${ANACONDA_ACTIVATE} ${PYTHON3_ENVIRONMENT}

python <<EOF
from tvb.recon.io.sensor import generate_schema_txt
generate_schema_txt("$1", "$PWD", "$2")
EOF