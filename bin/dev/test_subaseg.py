import os
(FREESURFER_HOME, SUBJECTS_DIR, SUBJECT, SUBJ_DIR, DMR, SEGMENT, ASEG_LIST,    
ASEG_LIST_lh, ASEG_LIST_rh, SURF_VN, VOL_VN, SURF, LABEL, MRI, SUBAPARC_AREA, SUBAPARC_MODE,  
T1_VOX2RASTKR_PATH, CRAS_PATH, VOX, STRMLNS_SIFT_NO) = [os.environ[key] for key in 'FREESURFER_HOME SUBJECTS_DIR \
    SUBJECT SUBJ_DIR DMR SEGMENT ASEG_LIST ASEG_LIST_lh ASEG_LIST_rh SURF_VN  \
    VOL_VN SURF LABEL MRI SUBAPARC_AREA SUBAPARC_MODE T1_VOX2RASTKR_PATH CRAS_PATH VOX STRMLNS_SIFT_NO'.split()]
#ASEG_LIST="8 10 11 12 13 16 17 18 26 47 49 50 51 52 53 54 58"
#ASEG_SURF="/Users/dionperd/CBR/VEP/JUNG/JUNG/aseg_surf"
import bnm.recon.algo.reconutils
import nibabel
import numpy as np


#for h in ['lh', 'rh']:
#    bnm.recon.algo.reconutils.sample_vol_on_surf(SURF+'/'+h+'.white',
#                                  SEGMENT+'/aparc+aseg-mask.nii.gz',
#                                  LABEL+'/'+h+'.aparc.annot',
#                                  SEGMENT+'/'+h+'.white-mask',
#                                  CRAS_PATH,ctx=h,
#                                  vn=int(SURF_VN),add_lbl=[2,41])
#
#for h in ['lh', 'rh']:
#    bnm.recon.algo.reconutils.sample_vol_on_surf(SURF+'/'+h+'.aseg',
#                                  SEGMENT+'/aparc+aseg-mask.nii.gz',
#                                  LABEL+'/'+h+'.aseg.annot',
#                                  SEGMENT+'/'+h+'.aseg-mask',
#                                  CRAS_PATH,ctx=None,
#                                  vn=int(SURF_VN),add_lbl=[])

SUBAPARC_MODE="con+geod+adj"
#SUBAPARC_MODE="geod"

ref_vol_path=SEGMENT+'/tdi_lbl-v'+VOX+'-in-t1.nii.gz'
consim_path=SEGMENT+'/consim-vol-counts'+STRMLNS_SIFT_NO+'-v'+VOX+'.npy'
labels=dict()
labels['lh']="1001 1007 1008 1009 1015 1030"
labels['rh']="2030"
#for h in ['lh', 'rh']:
h='rh'
bnm.recon.algo.reconutils.connectivity_geodesic_subparc(SURF+'/'+h+'.white',
                                             LABEL+'/'+h+'.aparc.annot',
                                             SEGMENT+'/'+h+'.white-mask-idx.npy', 
                                             out_annot_path=LABEL+'/'+h+'.aparc'+SUBAPARC_AREA+'-'+SUBAPARC_MODE+'.annot',
                                             parc_area=int(SUBAPARC_AREA), 
                                             labels=labels[h], hemi=None, ctx=h,  # labels=None, hemi=h,
                                             mode=SUBAPARC_MODE, 
                                             cras_path=CRAS_PATH, 
                                             ref_vol_path=ref_vol_path, 
                                             consim_path=consim_path, 
                                             lut_path=os.path.join(FREESURFER_HOME,'FreeSurferColorLUT.txt'))
#h='lh'
#bnm.recon.algo.reconutils.connectivity_geodesic_subparc(SURF+'/'+h+'.aseg',
#                                         LABEL+'/'+h+'.aseg.annot',
#                                         SEGMENT+'/'+h+'.aseg-mask-idx.npy', 
#                                         out_annot_path=LABEL+'/'+h+'.aseg'+SUBAPARC_AREA+'-'+SUBAPARC_MODE+'.annot',
#                                         parc_area=int(SUBAPARC_AREA), 
#                                         labels=ASEG_LIST_lh, hemi=None, ctx=None,
#                                         mode=SUBAPARC_MODE, 
#                                         cras_path=CRAS_PATH, 
#                                         ref_vol_path=ref_vol_path, 
#                                         consim_path=consim_path, 
#                                         lut_path=os.path.join(FREESURFER_HOME,'FreeSurferColorLUT.txt'))
#h='rh'
#bnm.recon.algo.reconutils.connectivity_geodesic_subparc(SURF+'/'+h+'.aseg',
#                                         LABEL+'/'+h+'.aseg.annot',
#                                         SEGMENT+'/'+h+'.aseg-mask-idx.npy', 
#                                         out_annot_path=LABEL+'/'+h+'.aseg'+SUBAPARC_AREA+'-'+SUBAPARC_MODE+'.annot',
#                                         parc_area=int(SUBAPARC_AREA), 
#                                         labels=ASEG_LIST_rh, hemi=None, ctx=None,
#                                         mode=SUBAPARC_MODE, 
#                                         mode=SUBAPARC_MODE, 
#                                         cras_path=CRAS_PATH, 
#                                         ref_vol_path=ref_vol_path, 
#                                         consim_path=consim_path, 
#                                         lut_path=os.path.join(FREESURFER_HOME,'FreeSurferColorLUT.txt'))