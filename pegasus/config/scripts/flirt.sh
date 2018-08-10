#!/usr/bin/env bash

export FSL_DIR
export PATH=${FSL_DIR}/bin:${PATH}
source ${FSL_DIR}/etc/fslconf/fsl.sh

if [ $# -eq 5 ]
then

echo flirt -in $1 -ref $2 -omat $3 -out $4 -dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -cost $5
flirt -in $1 -ref $2 -omat $3 -out $4 -dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -cost $5

else

echo flirt -in $1 -ref $2 -omat $3 -out $4 -dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -cost mutualinfo
flirt -in $1 -ref $2 -omat $3 -out $4 -dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -cost mutualinfo

fi


