#!/bin/bash

# bootstrap cmake because distros never have up to date
# TODO mv this to openmeeg which is only dependee

curl -O https://cmake.org/files/v3.4/cmake-3.4.3.tar.gz
tar xzf cmake-3.4.3.tar.gz
pushd cmake-3.4.3
./configure --prefix=/work/env
make -j6 # 2m50s
make install
popd
