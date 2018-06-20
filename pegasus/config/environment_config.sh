#!/usr/bin/env bash

#TODO place these in a better form ?

export PEGASUSHOME=/home/submitter
export PEGASUSSUBMIT=$${PEGASUSHOME}/submit
export PEGASUSSCRATCH=$${PEGASUSHOME}/scratch

export OS=LINUX

export FREESURFER_HOME=/opt/freesurfer-stable/freesurfer
export SUBJECTS_DIR=$${FREESURFER_HOME}/subjects
export FUNCTIONALS_DIR=$${FREESURFER_HOME}/sessions

export FSL_DIR=/usr/local/fsl

export MRTRIX_BIN=/opt/mrtrix3/bin
export MRTRIX_SCRIPTS=/opt/mrtrix3/scripts

export ANACONDA_ACTIVATE=/opt/conda/bin/activate
export PYTHON3_ENVIRONMENT=tvb_recon_python3_env

export SH_CUSTOM_FILES=/opt/tvb-recon/pegasus/config/scripts

export MNE_ROOT=/opt/MNE-2.7.4-3378-MacOSX-x86_64/bin