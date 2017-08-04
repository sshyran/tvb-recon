#!/usr/bin/env bash

export HOME=/Users/pipeline
export PATH=/WORK/MRtrix/mrtrix3/release/bin:/WORK/MRtrix/mrtrix3/scripts:${PATH}

dwi2response tournier $1 $2 -mask $3