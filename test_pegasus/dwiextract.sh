#!/usr/bin/env bash

export HOME=/Users/pipeline
export PATH=/WORK/MRtrix/mrtrix3/release/bin:$PATH

dwiextract -bzero $1 $2