# This setup should be done only once, before starting to use the bnm.recon.qc commands
sudo python setup.py develop

# The minimum environment variables needed for surface transformation
export SNAPSHOT_NUMBER=0
export FIGS=$FIGS

# This is needed for log messages
mkdir output
cd output
touch bnm.log
cd ..

# Surface transformation command calls

# This applies the transformation matrix from the surface.gii metadata
# python -m tvb.recon.qc.surface_transform surface_path.gii output_surface_path.gii

# This allows you to give a list of transformation matrices which will be all applied
# python -m tvb.recon.qc.surface_transform surface_path.gii output_surface_path.gii -matrix_paths matrix1.txt matrix2.txt

# To view the snapshots of a surface with an annotation applied, use this
# python -m tvb.recon.qc surf_annot surface_path.gii annotation_path