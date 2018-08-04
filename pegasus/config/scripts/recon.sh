#!/usr/bin/env bash

export HOME
export FREESURFER_HOME
export SUBJECTS_DIR
source ${FREESURFER_HOME}/FreeSurferEnv.sh

f=$PWD

# We assume that every case of rerunning should just resume recon-all without overwriting
# TODO: a proper management of recon-all
if [ -d "${SUBJECTS_DIR}/$1" ]; then
    echo OVERWRITE_RECONALL_FLAG=$5
    # if the previous recon-all run hasn't finished (i.e., generated wmparc.mgz) or overwrite flag is True, rerun/resume
    if [ ! -f "${SUBJECTS_DIR}/$1/mri/wmparc.mgz" ] || [ $5 == "True" ];
    then
        for h in "lh rh lh+rh"
        do
            if [ -f "${SUBJECTS_DIR}/$1/scripts/IsRunning.$h" ]; then
                rm ${SUBJECTS_DIR}/$1/scripts/IsRunning.$h
            fi
        done
        echo Restarting or resuming recon-all!
        recon-all -make all -no-isrunning -parallel -openmp $3 -s $1
    else
        # ...otherwise, skip recon-all!
        echo Skipping recon-all!
    fi
else
    recon-all -all -parallel -openmp $3 -s $1 -i $2
fi

cd ${SUBJECTS_DIR}/$1/mri
echo Copying T1.mgz, norm.mgz and brain.mgz
cp T1.mgz $f
cp norm.mgz $f
cp brain.mgz $f
cp aparc+aseg.mgz $f
for atlas_suffix in $4
do
    echo Copying aparc${atlas_suffix}+aseg.mgz
    cp aparc${atlas_suffix}+aseg.mgz $f
done

cd ../surf
echo Copying l/rh.pial and l/rh.white
cp lh.pial $f
cp rh.pial $f
cp lh.white $f
cp rh.white $f

cd ../label
cp lh.aparc.annot $f
cp rh.aparc.annot $f
for atlas_suffix in $4
do
    echo Copying l/rh.aparc${atlas_suffix}.annot
    cp lh.aparc${atlas_suffix}.annot $f
    cp rh.aparc${atlas_suffix}.annot $f
done