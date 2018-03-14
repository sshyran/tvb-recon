#!/bin/bash

# Workflow for generating forward models (gain matrices) for s/M/EEG.

# TODO visualization for BEM surfaces, sources & sensors
# TODO generate subcortical sources
# TODO average pial and white surface for cortical generator surface

# move to subject topic folder
pushd $BEM

for h in rh lh;
do
    # convert cortical surfaces format
    cp $SURF/$h.white.fsaverage5 ./cortical-$h
    python -c "import tvb.recon.algo.reconutils; tvb.recon.algo.reconutils.convert_fs_to_brain_visa('cortical-$h')"
    # source model for cortical hemispheres
    om_assemble -SurfSourceMat head_model.{geom,cond} cortical-$h.{tri,ssm}
done

# TODO: make source model for cortical dipoles
python -c "import tvb.recon.algo.reconutils; tvb.recon.algo.reconutils.gen_cort_sources()"
#Cortical volume sources model
for h in rh lh
do
    om_assemble -DipSourceMat head_model.{geom,cond} ./cortical-$h.{dip,dsm} # 2m32s
done

# TODO: make source model for subcortical
python -c "import tvb.recon.algo.reconutils; tvb.recon.algo.reconutils.gen_subcort_sources()"
#Subcortical volume sources model
for h in rh lh
do
    om_assemble -DipSourceMat head_model.{geom,cond} ./subcortical-$h.{dip,dsm} # 2m32s
done


popd