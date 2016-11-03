#!/bin/bash

DIR=$(cd $(dirname $0) && pwd)

if [ $# -ne 1 ]; then
    echo "Usage: $0 DAXFILE"
    exit 1
fi

DAXFILE=$1

PEGASUS_HOME=/opt/pegasus-home
export PEGASUS_HOME

OS=MACOSX
export OS

# This command tells Pegasus to plan the workflow contained in 
# dax file passed as an argument. The planned workflow will be stored
# in the "submit" directory. The execution site is "local".

pegasus-plan --conf pegasus.properties \
    --dax ${DAXFILE} \
    --dir ${PEGASUS_HOME}/submit \
    --cleanup leaf \
    --force \
    --sites condorpool \
    --submit
