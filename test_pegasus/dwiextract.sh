#!/usr/bin/env bash

export PATH="/WORK/MRtrix/mrtrix3/release/bin:/WORK/MRtrix/mrtrix3/scripts:$PATH"

/WORK/MRtrix/mrtrix3/release/bin/dwiextract -bzero $1 $2