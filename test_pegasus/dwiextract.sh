#!/usr/bin/env bash

export HOME
export PATH=${MRTRIX_BIN}:$PATH

dwiextract -bzero $1 $2