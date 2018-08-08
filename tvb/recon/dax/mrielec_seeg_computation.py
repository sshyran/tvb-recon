# -*- coding: utf-8 -*-

from Pegasus.DAX3 import File, Job, Link
from tvb.recon.dax.mappings import Inputs, MriElecSEEGCompFiles, MriElecSEEGCompJobNames, T1Files, CoregJobNames, \
    SEEGCompJobNames


class MriElecSEEGComputation(object):
    def __init__(self, subj, mri_elec_frmt, same_space_vol_pom, intensity_th, dilate=0, erode=0):
        self.subj = subj
        self.mri_elec_frmt = mri_elec_frmt
        self.same_space_vol_pom = same_space_vol_pom
        self.intensity_th = intensity_th
        self.dilate = dilate
        self.erode = erode

    def add_computation_steps(self, dax):
        mri_elec_nii_gz = File(Inputs.MRI_ELEC_INPUT.value)
        elecs_pom = File(Inputs.ELECS_POM.value)

        seeg_pom_xyz_txt = File(MriElecSEEGCompFiles.SEEG_POM_XYZ_TXT.value)
        seeg_pom_vol_nii_gz = File(MriElecSEEGCompFiles.SEEG_POM_VOL_NII_GZ.value)

        job1 = Job(MriElecSEEGCompJobNames.EXTRACT_POSITIONS_FROM_POM.value)
        job1.addArguments(elecs_pom, seeg_pom_xyz_txt, mri_elec_nii_gz, seeg_pom_vol_nii_gz)
        job1.uses(elecs_pom, Link.INPUT)
        job1.uses(mri_elec_nii_gz, Link.INPUT)
        job1.uses(seeg_pom_xyz_txt, Link.OUTPUT, transfer=True, register=True)
        job1.uses(seeg_pom_vol_nii_gz, Link.OUTPUT, transfer=True, register=True)

        dax.addJob(job1)

        t1_nii_gz = File(T1Files.T1_NII_GZ.value)
        mri_elec_to_t1_mat = File(MriElecSEEGCompFiles.MRI_ELEC_TO_T1_MAT.value)
        mri_elec_in_t1_nii_gz = File(MriElecSEEGCompFiles.MRI_ELEC_IN_T1_NII_GZ.value)
        job2 = Job(CoregJobNames.FLIRT.value)
        job2.addArguments(mri_elec_nii_gz, t1_nii_gz, mri_elec_to_t1_mat, mri_elec_in_t1_nii_gz)
        job2.uses(mri_elec_nii_gz, Link.INPUT)
        job2.uses(t1_nii_gz, Link.INPUT)
        job2.uses(mri_elec_to_t1_mat, Link.OUTPUT, transfer=True, register=True)
        job2.uses(mri_elec_in_t1_nii_gz, Link.OUTPUT, transfer=True, register=True)

        dax.addJob(job2)
        dax.depends(job2, job1)

        last_job = job2

        seeg_xyz_txt = File(MriElecSEEGCompFiles.SEEG_XYZ_TXT.value)
        seeg_in_t1_nii_gz = File(MriElecSEEGCompFiles.SEEG_IN_T1_NII_GZ.value)

        if self.same_space_vol_pom == "True":
            coords_pom_xyz_txt = File(MriElecSEEGCompFiles.COORDS_POM_XYZ_TXT.value)
            coords_xyz_txt = File(MriElecSEEGCompFiles.COORDS_XYZ_TXT.value)
            job3 = Job(MriElecSEEGCompJobNames.TRANSFORM_MERIELEC_COORDS.value)
            job3.addArguments(seeg_pom_xyz_txt, mri_elec_nii_gz, t1_nii_gz, seeg_xyz_txt, seeg_in_t1_nii_gz,
                              mri_elec_to_t1_mat, "")
            job3.uses(seeg_pom_xyz_txt, Link.INPUT)
            job3.uses(mri_elec_nii_gz, Link.INPUT)
            job3.uses(t1_nii_gz, Link.INPUT)
            job3.uses(mri_elec_to_t1_mat, Link.INPUT)
            job3.uses(seeg_xyz_txt, Link.OUTPUT, transfer=True, register=True)
            job3.uses(seeg_in_t1_nii_gz, Link.OUTPUT, transfer=True, register=True)
            job3.uses(coords_pom_xyz_txt, Link.OUTPUT, transfer=True, register=True)
            job3.uses(coords_xyz_txt, Link.OUTPUT, transfer=True, register=True)

            dax.addJob(job3)
            dax.depends(job3, job2)

            last_job = job3

        else:
            mri_elec_seeg = "mrielec_seeg.nii.gz"
            job3 = Job(SEEGCompJobNames.MRI_BINARIZE.value)

            #TODO: instead of 1 should be maxvoxel of mrielec volume
            arguments_list = []
            if self.dilate > 0:
                arguments_list.append("--dilate")
                arguments_list.append(str(self.dilate))
            if self.erode > 0:
                arguments_list.append("--erode")
                arguments_list.append(str(self.erode))

            job3.addArguments("--i", mri_elec_nii_gz, "--o", mri_elec_seeg, "--min", "1")
            job3.uses(mri_elec_nii_gz, Link.INPUT)
            job3.uses(mri_elec_seeg, Link.OUTPUT, transfer=True, register=True)

            dax.addJob(job3)
            dax.depends(job3, job2)


            seeg_pom_vol_bin = "seeg_pom_vol_bin.nii.gz"
            job4 = Job(SEEGCompJobNames.MRI_BINARIZE.value)

            arguments_list = []
            #TODO: same erode and dilate vals should be used for both steps?
            if self.dilate > 0:
                arguments_list.append("--dilate")
                arguments_list.append(str(self.dilate))
            if self.erode > 0:
                arguments_list.append("--erode")
                arguments_list.append(str(self.erode))

            job4.addArguments("--i", seeg_pom_vol_nii_gz, "--o", seeg_pom_vol_bin, "--min", "1")
            job4.uses(seeg_pom_vol_nii_gz, Link.INPUT)
            job4.uses(seeg_pom_vol_bin, Link.OUTPUT, transfer=True, register=True)

            dax.addJob(job4)
            dax.depends(job4, job3)


            pom_to_mri_elec = "pom_to_mrielec.mat"
            pom_in_mri_elec_seeg = "pom_in_mrielec_seeg.nii.gz"

            job5 = Job(CoregJobNames.FLIRT.value)
            job5.addArguments(seeg_pom_vol_bin, mri_elec_seeg, pom_to_mri_elec, pom_in_mri_elec_seeg)
            job5.uses(seeg_pom_vol_bin, Link.INPUT)
            job5.uses(mri_elec_seeg, Link.INPUT)
            job5.uses(pom_to_mri_elec, Link.OUTPUT, transfer=True, register=True)
            job5.uses(pom_in_mri_elec_seeg, Link.OUTPUT, transfer=True, register=True)

            dax.addJob(job5)
            dax.depends(job5, job4)


            seeg_pom_xyz_in_mrielec_xyz = MriElecSEEGCompFiles.SEEG_POM_XYZ_IN_MRIELEC_XYZ.value
            seeg_pom_vol_in_mrielec_nii_gz = MriElecSEEGCompFiles.SEEG_POM_VOL_IN_MRIELEC_NII_GZ.value

            job6 = Job(MriElecSEEGCompJobNames.TRANSFORM_MERIELEC_COORDS.value)
            job6.addArguments(seeg_pom_xyz_txt, seeg_pom_vol_nii_gz, mri_elec_nii_gz, seeg_pom_xyz_in_mrielec_xyz,
                              seeg_pom_vol_in_mrielec_nii_gz, pom_to_mri_elec)
            job6.uses(seeg_pom_xyz_txt, Link.INPUT)
            job6.uses(seeg_pom_vol_nii_gz, Link.INPUT)
            job6.uses(mri_elec_nii_gz, Link.INPUT)
            job6.uses(seeg_pom_xyz_in_mrielec_xyz, Link.OUTPUT, transfer=True, register=True)
            job6.uses(seeg_pom_vol_in_mrielec_nii_gz, Link.OUTPUT, transfer=True, register=True)
            job6.uses(pom_to_mri_elec, Link.INPUT)

            dax.addJob(job6)
            dax.depends(job6, job3)

            last_job = job6

            trsfrm_job = Job(MriElecSEEGCompJobNames.TRANSFORM_MERIELEC_COORDS.value)
            trsfrm_job.addArguments(seeg_pom_xyz_in_mrielec_xyz, seeg_pom_vol_in_mrielec_nii_gz, t1_nii_gz, seeg_xyz_txt,
                              seeg_in_t1_nii_gz, mri_elec_to_t1_mat)
            trsfrm_job.uses(seeg_pom_xyz_in_mrielec_xyz, Link.INPUT)
            trsfrm_job.uses(seeg_pom_vol_in_mrielec_nii_gz, Link.INPUT)
            trsfrm_job.uses(t1_nii_gz, Link.INPUT)
            trsfrm_job.uses(seeg_xyz_txt, Link.OUTPUT, transfer=True, register=True)
            trsfrm_job.uses(seeg_in_t1_nii_gz, Link.OUTPUT, transfer=True, register=True)
            trsfrm_job.uses(mri_elec_to_t1_mat, Link.INPUT)

            dax.addJob(trsfrm_job)
            dax.depends(trsfrm_job, last_job)
            last_job = trsfrm_job

        return last_job
