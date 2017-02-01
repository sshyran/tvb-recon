#!/usr/bin/env bash

set -eu
set -o pipefail

export PREFIX=${PREFIX:-"/work/env"}

if [[ -f $PREFIX/bin/cmake ]]
then
    echo "[30-cmake.sh] found our cmake, not building."
    exit 0
else
    echo "[30-cmake.sh] building newer cmake."
fi

# bootstrap cmake because distros never have up to date
curl -O https://cmake.org/files/v3.4/cmake-3.4.3.tar.gz
tar xzf cmake-3.4.3.tar.gz
pushd cmake-3.4.3
if [[ -z $(which cmake) ]]
then
    ./configure --prefix=$PREFIX
else
    cmake -DCMAKE_INSTALL_PREFIX=$PREFIX .
    make -j6
    make install
fi
popd
