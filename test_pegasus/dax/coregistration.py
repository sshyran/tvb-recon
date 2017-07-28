from Pegasus.DAX3 import File, Job, Link
from mappings import CoregFiles, CoregJobNames, T1Files, DWIFiles


class Coregistration(object):
    def __init__(self, use_flirt=True):
        self.use_flirt = use_flirt

    def _add_flirt_steps(self, dax, job_b0, job_t1, job_aparc_aseg):
        b0_nii_gz = File(DWIFiles.B0_NII_GZ.value)
        t1_nii_gz = File(T1Files.T1_NII_GZ.value)
        d2t_mat = File(CoregFiles.D2T_MAT.value)
        b0_in_t1 = File(CoregFiles.B0_IN_T1.value)
        job1 = Job(CoregJobNames.FLIRT.value, node_label="Register DWI to T1")
        job1.addArguments(b0_nii_gz, t1_nii_gz, d2t_mat, b0_in_t1)
        job1.uses(b0_nii_gz, link=Link.INPUT)
        job1.uses(t1_nii_gz, link=Link.INPUT)
        job1.uses(d2t_mat, link=Link.OUTPUT, transfer=False, register=False)
        job1.uses(b0_in_t1, link=Link.OUTPUT, transfer=True, register=False)
        dax.addJob(job1)

        dax.depends(job1, job_t1)
        dax.depends(job1, job_b0)

        t2d_mat = File(CoregFiles.T2D_MAT.value)
        job2 = Job(CoregJobNames.CONVERT_XFM.value, node_label="Convert d2t matrix to t2d matrix")
        job2.addArguments("-omat", t2d_mat, "-inverse", d2t_mat)
        job2.uses(d2t_mat, link=Link.INPUT)
        job2.uses(t2d_mat, link=Link.OUTPUT, transfer=True, register=False)
        dax.addJob(job2)

        dax.depends(job2, job1)

        t1_in_d_nii_gz = File(CoregFiles.T1_IN_D.value)
        job3 = Job(CoregJobNames.FLIRT_REVERSED.value, node_label="Register T1 to DWI")
        job3.addArguments(t1_nii_gz, b0_nii_gz, t1_in_d_nii_gz, t2d_mat)
        job3.uses(t1_nii_gz, link=Link.INPUT)
        job3.uses(b0_nii_gz, link=Link.INPUT)
        job3.uses(t2d_mat, link=Link.INPUT)
        job3.uses(t1_in_d_nii_gz, link=Link.OUTPUT, transfer=True, register=False)
        dax.addJob(job3)

        dax.depends(job3, job2)

        aparc_aseg_nii_gz = File(T1Files.APARC_ASEG_NII_GZ.value)
        aparc_aseg_in_d_nii_gz = File(CoregFiles.APARC_AGEG_IN_D.value)
        job4 = Job(CoregJobNames.FLIRT_REVERSED.value, node_label="Register APARC+ASEG to DWI")
        job4.addArguments(aparc_aseg_nii_gz, b0_nii_gz, aparc_aseg_in_d_nii_gz, t2d_mat)
        job4.uses(aparc_aseg_nii_gz, link=Link.INPUT)
        job4.uses(b0_nii_gz, link=Link.INPUT)
        job4.uses(t2d_mat, link=Link.INPUT)
        job4.uses(aparc_aseg_in_d_nii_gz, link=Link.OUTPUT, transfer=True, register=False)
        dax.addJob(job4)

        dax.depends(job4, job2)
        dax.depends(job4, job_aparc_aseg)

        return job3, job4

    def _add_fs_steps(self, dax, job_b0, job_t1, job_aparc_aseg):
        return None, None

    # job_b0 = job6, job_t1 = job7, job_aparc_aseg = job9
    def add_coregistration_steps(self, dax, job_b0, job_t1, job_aparc_aseg):

        if self.use_flirt is True:
            return self._add_flirt_steps(dax, job_b0, job_t1, job_aparc_aseg)
        else:
            return self._add_fs_steps(dax, job_b0, job_t1, job_aparc_aseg)
