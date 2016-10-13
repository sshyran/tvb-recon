#!/bin/bash

curl -L http://downloads.sourceforge.net/project/matio/matio/1.5.6/matio-1.5.6.tar.gz > matio-1.5.6.tar.gz
tar xzf matio-1.5.6.tar.gz
pushd matio-1.5.6
./configure
make
make install
ldconfig
