#!/bin/bash

set -e

pushd /usr/local
tar xzf ~/Downloads/MNE-*-Linux-x86_64.tar.gz
mv MNE-* mne
popd
