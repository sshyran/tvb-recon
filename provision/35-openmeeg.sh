#!/bin/bash

set -eu
set -o pipefail

export PREFIX=${PREFIX:-"/work/env"}

if [[ ! -z $(which om_assemble) ]]
then
    echo "[35-openmeeg.sh] found om_assemble, not building OpenMEEG"
    exit 0
else
    echo "[35-openmeeg.sh] building OpenMEEG."
fi

git clone https://github.com/openmeeg/openmeeg

pushd openmeeg

mkdir build && cd build
cmake -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release \
    -DENABLE_PYTHON=OFF -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DBUILD_DOCUMENTATION=OFF -DBUILD_TUTORIALS=OFF ..

make -j6
make install

popd