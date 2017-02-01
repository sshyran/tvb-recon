#!/usr/bin/env bash

set -eu
set -o pipefail

export PREFIX=${PREFIX:-"/work/env"}

if [[ ! -z $(which tckgen) ]]
then
    echo "[40-mrtrix.sh] found tckgen, not rebuilding."
    exit 0
else
    echo "[40-mrtrix.sh] building mrtrix3."
fi

if [[ ! -d mrtrix3 ]]
then
    git clone https://github.com/mrtrix3/mrtrix3
fi

pushd mrtrix3
    git checkout 0.3.15
    ./configure -nogui
    ./build
popd
