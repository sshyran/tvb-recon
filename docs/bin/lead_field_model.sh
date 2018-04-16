#!/bin/bash

pushd $BEM

for h in lh rh;
do
    #SEEG:
    #cortical dipoles
    om_gain -InternalPotential ./head-inv.mat ./cortical.dsm ./seeg.h2ipm ./seeg-cortical-$h.ds2ipm ./seeg-cortical-$h_gain.mat
    #subcortical dipoles
    om_gain -InternalPotential ./head-inv.mat ./subcortical.dsm ./seeg.h2ipm ./seeg-subcortical-$h.ds2ipm ./seeg-subcortical-$h_gain.mat

    #MEG:
    #cortical surfaces
    #om_gain -MEG ./head-inv.mat ./cortical.ssm ./meg.h2mm ./meg-$h.ss2mm ./meg-cortical-$h_gain.mat
    #subcortical dipoles
    #om_gain -MEG ./head-inv.mat ./subcortical.dsm ./meg.h2mm ./meg-$h.ds2mm ./meg-subcortical-$h_gain.mat

    #EEG:
    #cortical surfaces
    #om_gain -EEG ./head-inv.mat ./cortical.ssm ./EEG.h2em ./EEG-cortical-$h-gain.mat
    #subcortical dipoles
    #om_gain -EEG ./head-inv.mat ./subcortical.dsm ./EEG.h2em ./EEG-subcortical-$h-gain.mat
done


#TODO: merge/concatenate gain matrices to give seeg_gain.mat (similarly for EEG and MEG)

#Visual check:
# plot gain matrix
python<<EOF
import h5py, pylab as pl, numpy as np
linop = h5py.File('seeg_gain.mat')['/linop'][:]
print(linop.shape)
pl.hist(np.log(np.abs(linop.ravel())),100)
pl.show()
EOF

popd

#
#	    ctx-lh	ctx-rh	subcort-lh subcort-rh
#SEEG
#EEG
#MEG
#
#

