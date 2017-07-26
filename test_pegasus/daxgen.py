#!/usr/bin/env python
import os
import sys

import time
from Pegasus.DAX3 import ADAG, File, Job, Link, Executable, PFN, Profile, Namespace

if len(sys.argv) != 2:
    sys.stderr.write("Usage: %s DAXFILE\n" % (sys.argv[0]))
    sys.exit(1)
daxfile = sys.argv[1]

dax = ADAG("convert")
dax.metadata("created", time.ctime())

t1_input = File("t1_input.nii.gz")
t1_output = File("t1.mgz")
job1 = Job("mri_convert", node_label="T1 input conversion to MGZ")
job1.addArguments(t1_input, t1_output)
job1.uses(t1_input, link=Link.INPUT)
job1.uses(t1_output, link=Link.OUTPUT, transfer=True)
dax.addJob(job1)

# subject_dir = File("TVB2PEG2")
job2 = Job("recon", node_label="Recon-all for T1")
job2.addArguments("TVB2PEG22", t1_output)
job2.uses(t1_output, link=Link.INPUT)
# job2.uses(subject_dir, link=Link.OUTPUT, transfer=False, register=False)
dax.addJob(job2)

dax.depends(job2, job1)

t1_mgz_vol = File("T1.mgz")
t1_nii_gz_vol = File("T1.nii.gz")
job7 = Job("mri_convert", node_label="Convert T1 to NIFTI with good orientation")
job7.addArguments(t1_mgz_vol, t1_nii_gz_vol, "--out_orientation", "RAS")
job7.uses(t1_mgz_vol, link=Link.INPUT)
job7.uses(t1_nii_gz_vol, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job7)

dax.depends(job7, job2)

# wm_mgz_vol = File("wm.mgz")
# wm_nii_gz_vol = File("wm.nii.gz")
# job8 = Job("mri_convert", node_label="Convert WM to NIFTI with good orientation")
# job8.addArguments(wm_mgz_vol, wm_nii_gz_vol, "--out_orientation", "RAS")
# job8.uses(wm_mgz_vol, link=Link.INPUT)
# job8.uses(wm_nii_gz_vol, link=Link.OUTPUT, transfer=True, register=False)
# dax.addJob(job8)
#
# dax.depends(job8, job2)

aparc_aseg_mgz_vol = File("aparc+aseg.mgz")
aparc_aseg_nii_gz_vol = File("aparc+aseg.nii.gz")
job9 = Job("mri_convert", node_label="Convert APARC+ASEG to NIFTI with good orientation")
job9.addArguments(aparc_aseg_mgz_vol, aparc_aseg_nii_gz_vol, "--out_orientation", "RAS", "-rt", "nearest")
job9.uses(aparc_aseg_mgz_vol, link=Link.INPUT)
job9.uses(aparc_aseg_nii_gz_vol, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job9)

dax.depends(job9, job2)

# aseg_mgz_vol = File("aseg.mgz")
# aseg_nii_gz_vol = File("aseg.nii.gz")
# job22 = Job("mri_convert", node_label="Convert ASEG to NIFTI with good orientation")
# job22.addArguments(aseg_mgz_vol, aseg_nii_gz_vol, "--out_orientation", "RAS", "-rt", "nearest")
# job22.uses(aseg_mgz_vol, link=Link.INPUT)
# job22.uses(aseg_nii_gz_vol, link=Link.OUTPUT, transfer=True, register=False)
# dax.addJob(job22)
#
# dax.depends(job22, job2)

dwi_input = File("dwi_input.mif")
dwi_conv_output = File("dwi_raw.mif")
job3 = Job("mrconvert", node_label="Convert DWI to MIF")
job3.addArguments(dwi_input, dwi_conv_output, "-force")
job3.uses(dwi_input, link=Link.INPUT)
job3.uses(dwi_conv_output, link=Link.OUTPUT, transfer=False, register=False)
dax.addJob(job3)

dwi_pre_output = File("dwi.mif")
job4 = Job("dwipreproc", node_label="DWI preprocessing")
job4.addArguments("ap", dwi_conv_output, dwi_pre_output, "-rpe_none", "-nthreads", "2", "-force")
job4.uses(dwi_conv_output, link=Link.INPUT)
job4.uses(dwi_pre_output, link=Link.OUTPUT, transfer=False, register=False)
dax.addJob(job4)

dax.depends(job4, job3)

mask_output = File("mask.mif")
job5 = Job("dwi2mask", node_label="Create DWI mask")
job5.addArguments(dwi_pre_output, mask_output, "-nthreads", "2", "-force")
job5.uses(dwi_pre_output, link=Link.INPUT)
job5.uses(mask_output, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job5)

dax.depends(job5, job4)

b0_output = File("b0.nii.gz")
job6 = Job("dwiextract", node_label="Extract DWI B0")
job6.addArguments(dwi_pre_output, b0_output)
job6.uses(dwi_pre_output, link=Link.INPUT)
job6.uses(b0_output, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job6)

dax.depends(job6, job4)

d2t_mat = File("d2t.mat")
b0_in_t1 = File("b0-in-t1.nii.gz")
job10 = Job("flirt", node_label="Register DWI to T1")
job10.addArguments(b0_output, t1_nii_gz_vol, d2t_mat, b0_in_t1)
job10.uses(b0_output, link=Link.INPUT)
job10.uses(t1_nii_gz_vol, link=Link.INPUT)
job10.uses(d2t_mat, link=Link.OUTPUT, transfer=False, register=False)
job10.uses(b0_in_t1, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job10)

dax.depends(job10, job7)
dax.depends(job10, job6)

t2d_mat = File("t2d.mat")
job11 = Job("convert-xfm", node_label="Convert d2t matrix to t2d matrix")
job11.addArguments("-omat", t2d_mat, "-inverse", d2t_mat)
job11.uses(d2t_mat, link=Link.INPUT)
job11.uses(t2d_mat, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job11)

dax.depends(job11, job10)

t1_in_d_nii_gz = File("t1-in-d.nii.gz")
job12 = Job("flirt-reversed", node_label="Register T1 to DWI")
job12.addArguments(t1_nii_gz_vol, b0_output, t1_in_d_nii_gz, t2d_mat)
job12.uses(t1_nii_gz_vol, link=Link.INPUT)
job12.uses(b0_output, link=Link.INPUT)
job12.uses(t2d_mat, link=Link.INPUT)
job12.uses(t1_in_d_nii_gz, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job12)

dax.depends(job12, job11)

aparc_aseg_in_d_nii_gz = File("aparc+aseg-in-d.nii.gz")
job21 = Job("flirt-reversed", node_label="Register APARC+ASEG to DWI")
job21.addArguments(aparc_aseg_nii_gz_vol, b0_output, aparc_aseg_in_d_nii_gz, t2d_mat)
job21.uses(aparc_aseg_nii_gz_vol, link=Link.INPUT)
job21.uses(b0_output, link=Link.INPUT)
job21.uses(t2d_mat, link=Link.INPUT)
job21.uses(aparc_aseg_in_d_nii_gz, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job21)

dax.depends(job21, job9)

file_5tt = File("5tt.mif")
job13 = Job("ttgen", node_label="Generate 5tt MIF")
job13.addArguments(t1_in_d_nii_gz, file_5tt)
job13.uses(t1_in_d_nii_gz, link=Link.INPUT)
job13.uses(file_5tt, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job13)

dax.depends(job13, job12)

file_gmwmi = File("gmwmi.mif")
job14 = Job("tt2gmwmi", node_label="Extract GMWMI")
job14.addArguments(file_5tt, file_gmwmi, "-nthreads", "2")
job14.uses(file_5tt, link=Link.INPUT)
job14.uses(file_gmwmi, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job14)

dax.depends(job14, job13)

file_5ttvis = File("5ttvis.mif")
job15 = Job("tt2vis", node_label="Generate TT2VIS MIF")
job15.addArguments(file_5tt, file_5ttvis)
job15.uses(file_5tt, link=Link.INPUT)
job15.uses(file_5ttvis, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job15)

dax.depends(job15, job14)

file_respone = File("response.txt")
job16 = Job("dwi2response", node_label="Compute the DWI Response")
job16.addArguments(dwi_pre_output, file_respone, mask_output)
job16.uses(dwi_pre_output, link=Link.INPUT)
job16.uses(mask_output, link=Link.INPUT)
job16.uses(file_respone, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job16)

dax.depends(job16, job5)

file_wm_fod = File("wm_fod.mif")
job17 = Job("dwi2fod", node_label="Obtain WM FOD")
job17.addArguments("csd", dwi_pre_output, file_respone, file_wm_fod, "-mask", mask_output, "-nthreads", "2")
job17.uses(dwi_pre_output, link=Link.INPUT)
job17.uses(file_respone, link=Link.INPUT)
job17.uses(mask_output, link=Link.INPUT)
job17.uses(file_wm_fod, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job17)

dax.depends(job17, job16)

file_strmlns = File("25M.tck")
job18 = Job("tckgen", node_label="Generate tracts")
job18.addArguments(file_wm_fod, file_strmlns, "-number", "25M", "-seed_gmwmi", file_gmwmi, "-act", file_5tt,
                   "-unidirectional", "-maxlength", "250", "-step", "0.5", "-nthreads", "2")
job18.uses(file_wm_fod, link=Link.INPUT)
job18.uses(file_gmwmi, link=Link.INPUT)
job18.uses(file_5tt, link=Link.INPUT)
job18.uses(file_strmlns, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job18)

dax.depends(job18, job17)
dax.depends(job18, job13)

file_strmlns_sift = File("5M.tck")
job19 = Job("tcksift", node_label="Tracts SIFT")
job19.addArguments(file_strmlns, file_wm_fod, file_strmlns_sift, "-term_number", "5M", "-act", file_5tt, "-nthreads",
                   "2")
job19.uses(file_strmlns, link=Link.INPUT)
job19.uses(file_wm_fod, link=Link.INPUT)
job19.uses(file_5tt, link=Link.INPUT)
job19.uses(file_strmlns_sift, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job19)

dax.depends(job19, job18)
dax.depends(job19, job13)

file_tdi_ends = File("tdi_ends.mif")
job20 = Job("tckmap", node_label="TCKMAP")
job20.addArguments(file_strmlns_sift, file_tdi_ends, "-vox", "1", "-template", b0_output)
job20.uses(file_strmlns_sift, link=Link.INPUT)
job20.uses(b0_output, link=Link.INPUT)
job20.uses(file_tdi_ends, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job20)

dax.depends(job20, job19)

file_vol_lbl = File("aparc_aseg_lbl.nii.gz")
fs_color_lut = File("fs_color_lut.txt")
fs_default = File("fs_default.txt")
job23 = Job("labelconvert", node_label="Compute APARC+ASEG labeled for tracts")
job23.addArguments(aparc_aseg_in_d_nii_gz, fs_color_lut, fs_default, file_vol_lbl)
job23.uses(aparc_aseg_in_d_nii_gz, link=Link.INPUT)
job23.uses(fs_color_lut, link=Link.INPUT)
job23.uses(fs_default, link=Link.INPUT)
job23.uses(file_vol_lbl, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job23)

dax.depends(job23, job21)

file_aparc_aseg_counts5M_csv = File("aparc+aseg_counts5M.csv")
job24 = Job("tck2connectome", node_label="Generate weigths")
job24.addArguments(file_strmlns_sift, file_vol_lbl, "-assignment_radial_search", "2", file_aparc_aseg_counts5M_csv)
job24.uses(file_strmlns_sift, link=Link.INPUT)
job24.uses(file_vol_lbl, link=Link.INPUT)
job24.uses(file_aparc_aseg_counts5M_csv, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job24)

dax.depends(job24, job19)
dax.depends(job24, job23)

file_aparc_aseg_mean_tract_lengths5M_csv = File("aparc+aseg_mean_tract_lengths5M.csv")
job25 = Job("tck2connectome", node_label="Generate tract lengths")
job25.addArguments(file_strmlns_sift, file_vol_lbl, "-assignment_radial_search", "2", "-scale_length", "-stat_edge",
                   "mean", file_aparc_aseg_mean_tract_lengths5M_csv)
job25.uses(file_strmlns_sift, link=Link.INPUT)
job25.uses(file_vol_lbl, link=Link.INPUT)
job25.uses(file_aparc_aseg_mean_tract_lengths5M_csv, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job25)

dax.depends(job25, job19)
dax.depends(job25, job23)


f = open(daxfile, "w")
dax.writeXML(f)
f.close()
