#!/bin/bash

curl -O http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.16.tar.bz2
tar xjf hdf5-1.8.16.tar.bz2
pushd hdf5-1.8.16
./configure --prefix=/usr/local
make
make install
popd
