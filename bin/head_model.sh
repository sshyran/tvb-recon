#!/bin/bash

# Workflow for generating forward models (gain matrices) for s/M/EEG.

# TODO visualization for BEM surfaces, sources & sensors
# TODO generate subcortical sources
# TODO average pial and white surface for cortical generator surface

if [[ -z "$SUBJECT" ]]; then echo "SUBJECT?"; exit 1; fi

# move to subject topic folder
pushd $BEM

# orient sensors via parcellation, use default or config
# with a label->sensor map, compute rotation, then use ray intersect with
# scalp surface for EEG and inflate for MEG. Sensors should be equidistant
# to scalp surface while satisfying the label->sensor association.

# generate BEM surfaces
mne_watershed_bem

# decimate surfaces
for surf in ./watershed/*_surface
do
    mris_decimate -d 0.1 ${surf} ${surf}-low
done

#Visual check
python<<EOF
import os
subj = os.environ['SUBJECT']
TRGSUBJECT = os.environ['TRGSUBJECT']
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import nibabel.freesurfer
app = pg.mkQApp()
win = gl.GLViewWidget()
win.show()
def plotsurf(fname, **kwds):
    v, f = nibabel.freesurfer.read_geometry(fname)
    md = gl.MeshData(vertexes=v, faces=f)
    mi = gl.GLMeshItem(meshdata=md, drawEdges=True, drawFaces=False, **kwds)
    win.addItem(mi)
bem_templ = './watershed/{subj}_{surf}_surface-low'
for surf in 'brain inner_skull outer_skull outer_skin'.split():
    plotsurf(bem_templ.format(subj=subj, surf=surf))
plotsurf('../surf/lh.pial-'+TRGSUBJECT, color='r')
app.exec_()
EOF

# Visual check
#freeview -v $MRI/T1.mgz -f ./watershed/*_surface-low -viewport coronal -screenshot $FIGS/bem_surfs_low.png
python -m $SNAPSHOT --center_surface --snapshot_name watershed_lowsurf_t1 vol_surf $MRI/T1.nii.gz ./watershed/*_surface-low


# convert surfaces to BrainVisa format
for surf in ./watershed/*_surface-low
do
    python -c "import tvb.recon.algo.reconutils; tvb.recon.algo.reconutils.convert_fs_to_brain_visa('$surf')"
done

# build head matrix
python -c "import tvb.recon.algo.reconutils; tvb.recon.algo.reconutils.gen_head_model()"

om_assemble -HM ./head_model.geom ./head_model.cond ./head.mat # 2m32s
om_minverser ./head.mat ./head-inv.mat # 3m30s

popd
