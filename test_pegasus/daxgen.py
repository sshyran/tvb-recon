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

# t1_input = File("t1_input.nii.gz")
# t1_output = File("t1.mgz")
# job1 = Job("mri_convert")
# job1.addArguments(t1_input, t1_output)
# job1.uses(t1_input, link=Link.INPUT)
# job1.uses(t1_output, link=Link.OUTPUT, transfer=True)
# dax.addJob(job1)
#
# job2 = Job("recon")
# job2.addArguments("TVB2PEG22", t1_output)
# job2.uses(t1_output, link=Link.INPUT)
# dax.addJob(job2)
#
# dax.depends(job2, job1)

t1_mgz_vol = File("T1.mgz")
t1_nii_gz_vol = File("T1.nii.gz")
job7 = Job("mri_convert")
job7.addArguments(t1_mgz_vol, t1_nii_gz_vol, "--out_orientation", "RAS")
job7.uses(t1_mgz_vol, link=Link.INPUT)
job7.uses(t1_nii_gz_vol, link=Link.OUTPUT, transfer=True, register=False)

dax.addJob(job7)

wm_mgz_vol = File("wm.mgz")
wm_nii_gz_vol = File("wm.nii.gz")
job8 = Job("mri_convert")
job8.addArguments(wm_mgz_vol, wm_nii_gz_vol, "--out_orientation", "RAS")
job8.uses(wm_mgz_vol, link=Link.INPUT)
job8.uses(wm_nii_gz_vol, link=Link.OUTPUT, transfer=True, register=False)

dax.addJob(job8)

for aseg_file in "aparc+aseg", "aseg":
    mgz_vol = File(aseg_file + ".mgz")
    nii_gz_vol = File(aseg_file + ".nii.gz")
    job9 = Job("mri_convert")
    job9.addArguments(mgz_vol, nii_gz_vol, "--out_orientation", "RAS", "-rt", "nearest")
    job9.uses(mgz_vol, link=Link.INPUT)
    job9.uses(nii_gz_vol, link=Link.OUTPUT, transfer=True, register=False)

    dax.addJob(job9)

dwi_input = File("dwi_input.mif")
dwi_conv_output = File("dwi_raw.mif")
job3 = Job("mrconvert")
job3.addArguments(dwi_input, dwi_conv_output)
job3.uses(dwi_input, link=Link.INPUT)
job3.uses(dwi_conv_output, link=Link.OUTPUT, transfer=False, register=False)
dax.addJob(job3)

dwi_pre_output = File("dwi.mif")
job4 = Job("dwipreproc")
job4.addArguments("ap", dwi_conv_output, dwi_pre_output, "-rpe_none", "-nthreads", "2", "-force")
job4.uses(dwi_conv_output, link=Link.INPUT)
job4.uses(dwi_pre_output, link=Link.OUTPUT, transfer=False, register=False)
dax.addJob(job4)

dax.depends(job4, job3)

mask_output = File("mask.mif")
job5 = Job("dwi2mask")
job5.addArguments(dwi_pre_output, mask_output, "-nthreads", "2", "-force")
job5.uses(dwi_pre_output, link=Link.INPUT)
job5.uses(mask_output, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job5)

dax.depends(job5, job4)

b0_output = File("b0.nii.gz")
job6 = Job("dwiextract")
job6.addArguments(dwi_pre_output, b0_output)
job6.uses(dwi_pre_output, link=Link.INPUT)
job6.uses(b0_output, link=Link.OUTPUT, transfer=True, register=False)
dax.addJob(job6)

dax.depends(job6, job5)

d2t_mat = File("d2t.mat")
b0_in_t1 = File("b0-in-t1.nii.gz")
job10 = Job("flirt")
job10.addArguments("-in", b0_output, "-ref", t1_nii_gz_vol, "-omat", d2t_mat, "-out", b0_in_t1, "-dof", "12",
                   "-searchrx", "-180 180", "-searchry", "-180 180", "-searchrz", "-180 180", "-cost", "mutualinfo")
job10.uses(b0_output, link=Link.INPUT)
job10.uses(t1_nii_gz_vol, link=Link.INPUT)
job10.uses(d2t_mat, link=Link.OUTPUT, transfer=True, register=False)
job10.uses(b0_in_t1, link=Link.OUTPUT, transfer=False, register=False)
dax.addJob(job10)

dax.depends(job10, job7)
dax.depends(job10, job6)

f = open(daxfile, "w")
dax.writeXML(f)
f.close()
