#!/usr/bin/env bash

# not included via Vagrant's shell provisioner..
source /etc/bash.bashrc

set -eu
set -o pipefail

export PREFIX=${PREFIX:-"/work/env"}

if [[ -f $PREFIX/jupyter-lab.txt ]]
then
    echo "[10-work-env.sh] $PREFIX/jupyter-lab.txt exists, skipping work env setup"
    exit 0
fi

mkdir -p $PREFIX/src
pushd $PREFIX/src

    bash /vagrant/provision/30-cmake.sh
    bash /vagrant/provision/35-openmeeg.sh
    bash /vagrant/provision/40-mrtrix.sh
    bash /vagrant/provision/70-python.sh

popd # $PREFIX/src

# unpack freesurfer
if [[ -z $(which mri_convert) ]]
then
    echo "[10-work-env.sh] unpacking freesurfer.."
    tar -C $PREFIX -xzf /vagrant/provision/freesurfer*.tar.gz
fi
cp /vagrant/provision/freesurfer-license.txt $PREFIX/freesurfer/license.txt

# unpack MNE
if [[ -z $(which mne_watershed_bem) ]]
then
    echo "[10-work-env.sh] unpacking MNE.."
    tar -C $PREFIX -xzf /vagrant/provision/MNE*.tar.gz
    ln -s $PREFIX/MNE* $PREFIX/mne
fi

# start jupyter lab and copy out connection info
sudo systemctl start jupyter-lab
sleep 5
sudo journalctl -u jupyter-lab > $PREFIX/jupyter-lab.txt

cat $PREFIX/jupyter-lab.txt
