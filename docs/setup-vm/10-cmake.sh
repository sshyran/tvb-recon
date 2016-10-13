#!/bin/bash

curl -O https://cmake.org/files/v3.4/cmake-3.4.3.tar.gz
tar xzf cmake-3.4.3.tar.gz
pushd cmake-3.4.3
mkdir build && cd build
cmake ..  # 1m30s
make # 2m50s
make install
popd
