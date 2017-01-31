#!/bin/bash

# TODO this builds hdf5 matio, so maybe drop our builds of those
git clone https://github.com/maedoc/openmeeg

pushd openmeeg
git checkout zlib-1.2.11

mkdir build && cd build
cmake -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release \
    -DENABLE_PYTHON=OFF -DCMAKE_INSTALL_PREFIX=/work/env \
    -DBUILD_DOCUMENTATION=OFF -DBUILD_TUTORIALS=OFF ..

make -j6
make install

popd