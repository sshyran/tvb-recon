#!/bin/bash

# Workflow for generating forward models (gain matrices) for s/M/EEG.

# TODO move to a Python script.. filenames don't scale well to handle all
# possibilities. several pieces of data need to be labeled by parcellation,
# modality, etc..

# TODO check inputs
if [[ -z "$SUBJECT" ]]; then echo "SUBJECT?"; exit 1; fi

# orient sensors via parcellation, use default or config
# with a label->sensor map, compute rotation, then use ray intersect with
# scalp surface for EEG and inflate for MEG. Sensors should be equidistant
# to scalp surface while satisfying the label->sensor association.

# generate BEM surfaces
mne_watershed_bem

# decimate surfaces
for surf in subjects/tvb/bem/watershed/*_surface
do
    time mris_decimate -d 0.1 ${surf} ${surf}-low
done

# inspect
freeview -v subjects/tvb/mri/T1.mgz -f subjects/tvb/bem/watershed/*-low -viewport coronal

# convert them to BrainVisa format
for surf in *_surface-low
do
    python<<EOF
import utils
utils.convert_fs_to_brain_visa("$surf")
EOF
done

# build head matrix
python -c "import utils; utils.gen_head_model()"
#pushd ${SUBJECTS_DIR}/${SUBJECT}/bem
om_assemble -HM head_model.geom head_model.cond head.mat # 2m32s
om_minverser head.mat head-inv.mat # 3m30s
#popd

# convert cortical surfaces format
for h in rh lh; do
    cp ../surf/$h.pial.fsaverage5 ./cortical-$h
    python -c "import utils; utils.convert_fs_to_brain_visa('cortical-$h')"
done

# make source models
python -c "import utils; utils.gen_subcort_sources()"
om_assemble -DSM head_model.{geom,cond} $subcortical.{dip,dsm}
for h in rh lh; do
    om_assemble -SSM head_model.{geom,cond} cortical-$h.{tri,ssm}
done

# make sensor models
om_assemble -h2em head_model.{geom,cond} EEG.sensors EEG.h2em
om_assemble -h2mm head_model.{geom,cond} MEG.sensors MEG.h2mm

# make gain matrices per source / sensor pair
for source_model in subcortical.dsm cortical-{lh,rh}.ssm
do
    for sensor_mode in EEG MEG
    do
        om_gain -$sensor_mode head-inv.mat $source_model $sensor_mode.* \
            $sensor_mode.gain.mat
    done
done

# apply parcs to gains for per-parc forwards
for parc in $parcnames
do
    for sensors in EEG MEG
    do
        # TODO
        python -c "import utils; utils.parc_gain('$parc', '$sensors')"
    done
done
