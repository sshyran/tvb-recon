#!/bin/bash

for script in 20-gsl 25-mrtrix 30-hdf5 32-matio 35-openmeeg
do
	bash ${TVBVIRT}/setup/${script}.sh
done
