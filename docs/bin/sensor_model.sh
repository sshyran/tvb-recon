#!/bin/bash

# Workflow for generating forward models (gain matrices) for s/M/EEG.

# TODO visualization for BEM surfaces, sources & sensors
# TODO generate subcortical sources
# TODO average pial and white surface for cortical generator surface

# move to subject topic folder
pushd $BEM

seeg_sensor_file=$SEEG_XYZ_FILE

#Sensor models...

#head to...
#om_assemble -h2em ./head_model.{geom, cond} ./EEG_sensors ./EEG.h2em # 2m32s
#om_assemble -h2mm ./head_model.{geom, cond} ./MEG_sensors ./MEG.h2mm # 2m32s
#...and for SEEG:
#head to internal potential model
om_assemble -h2ipm head_model.{geom,cond} $seeg_sensor_file seeg.h2ipm

for h in lh rh;
do
    for source_model in cortical subcortical
    do
        #(sub)cortical dipoles sources to internal potential model
        om_assemble -ds2ipm ./head_model.{geom,cond} ./$source_model-$h.dsm $seeg_sensor_file ./seeg-$source_model-$h.ds2ipm
    done
    #subcortical dipoles sources to MEG
    #om_assemble -ds2mm ./head_model.{geom,cond} ./subcortical-$h.dsm $meg_sensor_file ./meg-$h.ds2mm
    #cortical surface sources to MEG
    #om_assemble -ss2mm ./head_model.{geom,cond} ./cortical-$h.ssm $meg_sensor_file ./meg-$h.ss2mm
done


popd