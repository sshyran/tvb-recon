#!/bin/bash

DIR=$(cd $(dirname $0) && pwd)

if [ $# -ne 3 ]; then
    echo "Usage: $0 DAXFILE PEGASUS_PROPS ENV_CONFIG"
    exit 1
fi

DAXFILE=$1
PEGASUS_PROPS=$2
ENV_CONFIG=$3

source ${ENV_CONFIG}

# This command tells Pegasus to plan the workflow contained in 
# dax file passed as an argument. The planned workflow will be stored
# in the "submit" directory. The execution site is "local".

pegasus-plan --conf ${PEGASUS_PROPS} \
    --dax ${DAXFILE} \
    --dir ${PEGASUSSUBMIT} \
    --output-site local \
    --sites condorpool \
    --submit \
    --cleanup leaf