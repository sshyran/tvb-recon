#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR/..

# Generate with Pegasus pythonA PI and patient configuration a DAX file
bash generate_dax.sh extra/patient1.properties

# Generate Image Graph (DOT file) from DAX file
pegasus-graphviz -o output/main_bnm.dot output/main_bnm.dax

# Convert DOT image into PNG for compatibility
dot output/main_bnm.dot -Tpng -o output/main_bnm.png

# Plan DAX execution
bash plan_dax.sh output/main_bnm.dax
