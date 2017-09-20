#!/usr/bin/env bash

#TODO place these in a better form ?

export PEGASUSHOME=/WORK/pegasus-4.7.4
export PEGASUSSUBMIT=${PEGASUSHOME}/submit
export PEGASUSSCRATCH=${PEGASUSHOME}/scratch

export OS=MACOSX

export HOME=/Users/pipeline

export FREESURFER_HOME=/WORK/FS_NEW/freesurfer
export SUBJECTS_DIR=${FREESURFER_HOME}/subjects
export FUNCTIONALS_DIR=${FREESURFER_HOME}/sessions

export FSL_DIR=/WORK/FSL/fsl

export MRTRIX_BIN=/WORK/MRtrix/mrtrix3/release/bin
export MRTRIX_SCRIPTS=/WORK/MRtrix/mrtrix3/scripts

export ANACONDA_ACTIVATE=//anaconda/bin/activate
export PYTHON3_ENVIRONMENT=tvb_recon_python3_env

export SH_CUSTOM_FILES=/WORK/BNM/tvb-recon/tvb-recon/pegasus/config/scripts

export MNE_ROOT=/WORK/MNE/MNE-2.7.0-3106-MacOSX-i386