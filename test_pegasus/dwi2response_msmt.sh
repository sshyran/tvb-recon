#!/usr/bin/env bash

export HOME=/Users/pipeline
export PATH=/WORK/MRtrix/mrtrix3/release/bin:/WORK/MRtrix/mrtrix3/scripts:${PATH}

dwi2response msmt_5tt $1 $2 $3 $4 $5 -voxels $6 -nthreads $7