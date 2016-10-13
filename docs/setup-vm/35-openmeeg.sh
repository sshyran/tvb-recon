#!/bin/bash

git clone https://github.com/openmeeg/openmeeg
pushd openmeeg && mkdir build && cd build
cmake -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release ..
make
make install
ldconfig
popd
