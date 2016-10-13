#!/bin/bash

curl -O ftp://ftp.gnu.org/gnu/gsl/gsl-2.1.tar.gz
tar xzf gsl-2.1.tar.gz && pushd gsl-2.1
./configure && make && make install
ldconfig
popd

