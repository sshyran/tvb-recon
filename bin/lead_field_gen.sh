#!/bin/bash

pushd $BEM

#
#	ctx-lh	ctx-rh	subcort-lh subcort-rh
#SEEG
#EEG
#MEG
#
#

# Parcellations to be used should be given as inputs.
#VOLS=$*
#Sensors as well
#sensors=$*

#Apply vols to gains for per-vol lead fields
#TODO
for vol in $VOLS
do
    for sensors in SEEG EEG MEG
    do
        #TODO
        python -c "import bnm.recon.algo.reconutils; bnm.recon.algo.reconutils.parc_gain('$MRI/$vol','./$sensors')"
    done
done

popd