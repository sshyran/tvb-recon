#!/usr/bin/env bash

export FREESURFER_HOME="/WORK/FS_NEW/freesurfer"
source "/WORK/FS_NEW/freesurfer/SetUpFreeSurfer.sh"

/WORK/FS_NEW/freesurfer/bin/recon-all -all -s $1 -i $2