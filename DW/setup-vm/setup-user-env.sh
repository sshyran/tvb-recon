#!/bin/bash

mkdir $HOME/subjects

cat >> ~/.bashrc <<EOF
export SUBJECTS_DIR=$HOME/subjects
export FREESURFER_HOME=/usr/local/freesurfer
export FSLDIR=/usr/local/fsl
export MNE_ROOT=/usr/local/mne
export PATH=${FSLDIR}/bin:/usr/local/mrtrix3/bin:/usr/local/mrtrix3/scripts:$PATH
source ${FREESURFER_HOME}/FreeSurferEnv.sh
source ${MNE_ROOT}/bin/mne_setup_sh
EOF
