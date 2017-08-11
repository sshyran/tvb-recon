#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR/..

pwd

# Generate with Pegasus python API and patient configuration a DAX file
bash pegasus/generate_dax.sh pegasus/dax_output/main_bnm.dax pegasus/config/patient_flow.properties

# Generate Image Graph (DOT file) from DAX file
pegasus-graphviz -o pegasus/dax_output/main_bnm.dot pegasus/dax_output/main_bnm.dax
 
# Convert DOT image into PNG for compatibility
dot pegasus/dax_output/main_bnm.dot -Tpng -o pegasus/dax_output/main_bnm.png
 
# Plan DAX execution
bash pegasus/plan_dax.sh pegasus/dax_output/main_bnm.dax pegasus/config/pegasus.properties pegasus/config/environment_config.sh