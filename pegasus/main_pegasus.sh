#!/usr/bin/env bash
# This file expects 2 path arguments: 
#   1. $1 = Path to the patient configurations folder (contains: environment_config.sh, patient_flow.properties, 
#                                                       pegasus.properties, rc_out.txt, rc.txt, sites.xml, tc.txt)
#   2. $2 = Path to the patient dax folder (contains: main_bnm.dax, main_bnm.dot, main_bnm.png)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR/..

pwd

# Generate with Pegasus python API and patient configuration a DAX file
bash pegasus/generate_dax.sh $2/main_bnm.dax $1/patient_flow.properties

# Generate Image Graph (DOT file) from DAX file
pegasus-graphviz -o $2/main_bnm.dot $2/main_bnm.dax
 
# Convert DOT image into PNG for compatibility
dot $2/main_bnm.dot -Tpng -o $2/main_bnm.png
 
# Plan DAX execution
bash pegasus/plan_dax.sh $2/main_bnm.dax $1/pegasus.properties $1/environment_config.sh
