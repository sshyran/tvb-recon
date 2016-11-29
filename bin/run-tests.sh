#!/bin/bash

# ignore extern and scripts under dev
ignores=""
for ign in extern bin/dev; do ignores="--ignore=$ign $ignores"; done

# TODO adapt for different datasets, etc
export SUBJECTS_DIR=`pwd`
export SUBJECT='bnm'

# maybe do coverage
cov=""
if [[ $COV == "yes" ]]; then cov="--cov=bnm"; fi

# run 'em
py.test --cov-config .coveragec $cov $ignores
