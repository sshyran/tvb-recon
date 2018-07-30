# -*- coding: utf-8 -*-

from Pegasus.DAX3 import File, Job, Link
from tvb.recon.dax.mappings import Inputs, MriElecSEEGCompFiles, MriElecSEEGCompJobNames, T1Files, CoregJobNames


class MriElecSEEGComputation(object):
    def __init__(self, subj, mri_elec_frmt, same_space_vol_pom, intensity_th):
        self.subj = subj
        self.mri_elec_frmt = mri_elec_frmt
        self.same_space_vol_pom = same_space_vol_pom
        self.intensity_th = intensity_th

    def add_computation_steps(self, dax):
        mri_elec_nii_gz = File(Inputs.MRI_ELEC_INPUT.value)
        elecs_pom = File(Inputs.ELECS_POM.value)

        seeg_pom_xyz_txt = File(MriElecSEEGCompFiles.SEEG_POM_XYZ_TXT.value)
        seeg_pom_vol_txt = File(MriElecSEEGCompFiles.SEEG_POM_VOL_NII_GZ.value)

        job1 = Job(MriElecSEEGCompJobNames.EXTRACT_POSITIONS_FROM_POM.value)
        job1.addArguments(elecs_pom, seeg_pom_xyz_txt, mri_elec_nii_gz, seeg_pom_vol_txt)
        job1.uses(elecs_pom, Link.INPUT)
        job1.uses(mri_elec_nii_gz, Link.INPUT)
        job1.uses(seeg_pom_xyz_txt, Link.OUTPUT, transfer=True, register=True)
        job1.uses(seeg_pom_vol_txt, Link.OUTPUT, transfer=True, register=True)

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
            seeg_pom_xyz_in_mrielec_xyz = File(MriElecSEEGCompFiles.SEEG_POM_XYZ_IN_MRIELEC_XYZ.value)
            seeg_pom_vol_in_mrielec_nii_gz = File(MriElecSEEGCompFiles.POM_VOL_IN_MRIELEC_NII_GZ.value)
            job3 = Job(MriElecSEEGCompJobNames.COREGISTER_ELEC_POM_AND_MRI.value)
            job3.addArguments(seeg_pom_vol_txt, mri_elec_nii_gz, seeg_pom_xyz_txt, "10", "2")
            job3.uses(seeg_pom_vol_txt, Link.INPUT)
            job3.uses(mri_elec_nii_gz, Link.INPUT)
            job3.uses(seeg_pom_xyz_txt, Link.INPUT)
            job3.uses(seeg_pom_xyz_in_mrielec_xyz, Link.OUTPUT, transfer=True, register=True)
            job3.uses(seeg_pom_vol_in_mrielec_nii_gz, Link.OUTPUT, transfer=True, register=True)

            dax.addJob(job3)
            dax.depends(job3, job2)

            coords_pom_xyz_in_mrielec_xyz = File(MriElecSEEGCompFiles.COORDS_POM_XYZ_IN_MRIELEC_XYZ.value)
            coords_pom_vol_in_mrielec_nii_gz = File(MriElecSEEGCompFiles.COORDS_POM_VOL_IN_MRIELEC_NII_GZ.value)

            job4 = Job(MriElecSEEGCompJobNames.TRANSFORM_MERIELEC_COORDS.value)
            job4.addArguments(seeg_pom_xyz_in_mrielec_xyz, seeg_pom_vol_in_mrielec_nii_gz, t1_nii_gz, seeg_xyz_txt,
                              seeg_in_t1_nii_gz, mri_elec_to_t1_mat)
            job4.uses(seeg_pom_xyz_in_mrielec_xyz, Link.INPUT)
            job4.uses(seeg_pom_vol_in_mrielec_nii_gz, Link.INPUT)
            job4.uses(t1_nii_gz, Link.INPUT)
            job4.uses(mri_elec_to_t1_mat, Link.INPUT)
            job4.uses(seeg_xyz_txt, Link.OUTPUT, transfer=True, register=True)
            job4.uses(seeg_in_t1_nii_gz, Link.OUTPUT, transfer=True, register=True)
            # job4.uses(coords_pom_xyz_in_mrielec_xyz, Link.OUTPUT, transfer=True, register=True)
            # job4.uses(coords_pom_vol_in_mrielec_nii_gz, Link.OUTPUT, transfer=True, register=True)

            dax.addJob(job4)
            dax.depends(job4, job3)

            last_job = job4

        return last_job
