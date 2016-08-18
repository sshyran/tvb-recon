#!/bin/bash

set -e

ark=freesurfer-Linux-centos6_x86_64-dev.tar.gz

echo "Downloading FreeSurfer"
curl -O ftp://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/dev/$ark

here=$(pwd)
pushd /usr/local/
echo "unpacking FreeSurfer in /usr/local"
tar xzf $here/$ark
popd

rm $ark
