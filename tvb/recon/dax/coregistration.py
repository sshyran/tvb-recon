from Pegasus.DAX3 import File, Job, Link
from tvb.recon.dax import AtlasSuffix
from tvb.recon.dax.mappings import DWIFiles, T1Files, CoregFiles, CoregJobNames, T1JobNames
from tvb.recon.dax.qc_snapshots import QCSnapshots


class Coregistration(object):
    def __init__(self, subject, use_flirt=True, atlas_suffixes=[AtlasSuffix.DEFAULT]):
        self.subject = subject
        self.use_flirt = use_flirt
        self.qc_snapshots = QCSnapshots.get_instance()
        self.atlas_suffixes = atlas_suffixes

    def _add_flirt_steps(self, dax, job_b0, job_t1, jobs_aparc_aseg):
        b0_nii_gz = File(DWIFiles.B0_NII_GZ.value)
        t1_nii_gz = File(T1Files.T1_NII_GZ.value)
        d2t_mat = File(CoregFiles.D2T_MAT.value)
        b0_in_t1 = File(CoregFiles.B0_IN_T1.value)
        job1 = Job(CoregJobNames.FLIRT.value, node_label="Register DWI to T1")
        job1.addArguments(b0_nii_gz, t1_nii_gz, d2t_mat, b0_in_t1)
        job1.uses(b0_nii_gz, link=Link.INPUT)
        job1.uses(t1_nii_gz, link=Link.INPUT)
        job1.uses(d2t_mat, link=Link.OUTPUT, transfer=True, register=True)
        job1.uses(b0_in_t1, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job1)

        dax.depends(job1, job_t1)
        dax.depends(job1, job_b0)

        self.qc_snapshots.add_2vols_snapshot_step(dax, [job1], t1_nii_gz, b0_in_t1, "b0_in_t1")

        t2d_mat = File(CoregFiles.T2D_MAT.value)
        job2 = Job(CoregJobNames.CONVERT_XFM.value, node_label="Convert d2t matrix to t2d matrix")
        job2.addArguments("-omat", t2d_mat, "-inverse", d2t_mat)
        job2.uses(d2t_mat, link=Link.INPUT)
        job2.uses(t2d_mat, link=Link.OUTPUT, transfer=False, register=False)
        dax.addJob(job2)

        dax.depends(job2, job1)

        t1_in_d_nii_gz = File(CoregFiles.T1_IN_D.value)
        job3 = Job(CoregJobNames.FLIRT_REVERSED.value, node_label="Register T1 to DWI")
        job3.addArguments(t1_nii_gz, b0_nii_gz, t1_in_d_nii_gz, t2d_mat)
        job3.uses(t1_nii_gz, link=Link.INPUT)
        job3.uses(b0_nii_gz, link=Link.INPUT)
        job3.uses(t2d_mat, link=Link.INPUT)
        job3.uses(t1_in_d_nii_gz, link=Link.OUTPUT, transfer=False, register=False)
        dax.addJob(job3)

        dax.depends(job3, job2)

        self.qc_snapshots.add_2vols_snapshot_step(dax, [job3], t1_in_d_nii_gz, b0_nii_gz, "t1_in_b0")

        jobs4 = []
        aparc_aseg_nii_gz = []
        aparc_aseg_in_d_nii_gz = []
        for atlas_suffix, job_aparc_aseg in zip(self.atlas_suffixes, jobs_aparc_aseg):
            aparc_aseg_nii_gz.append(File(T1Files.APARC_ASEG_NII_GZ.value % atlas_suffix))
            aparc_aseg_in_d_nii_gz.append(File(CoregFiles.APARC_ASEG_IN_D.value % atlas_suffix))
            jobs4.append(Job(CoregJobNames.FLIRT_REVERSED.value,
                            node_label="Register APARC%s+ASEG to DWI" % atlas_suffix))
            jobs4[-1].addArguments(aparc_aseg_nii_gz[-1], b0_nii_gz, aparc_aseg_in_d_nii_gz[-1], t2d_mat)
            jobs4[-1].uses(aparc_aseg_nii_gz[-1], link=Link.INPUT)
            jobs4[-1].uses(b0_nii_gz, link=Link.INPUT)
            jobs4[-1].uses(t2d_mat, link=Link.INPUT)
            jobs4[-1].uses(aparc_aseg_in_d_nii_gz[-1], link=Link.OUTPUT, transfer=False, register=False)
            dax.addJob(jobs4[-1])

            dax.depends(jobs4[-1], job2)
            dax.depends(jobs4[-1], job_aparc_aseg)

            self.qc_snapshots.add_2vols_snapshot_step(dax, [jobs4[-1]],
                                                      aparc_aseg_in_d_nii_gz[-1], b0_nii_gz,
                                                      "aparc%s_aseg_in_b0" % atlas_suffix)
            self.qc_snapshots.add_3vols_snapshot_step(dax, [job3, jobs4[-1]], t1_in_d_nii_gz, b0_nii_gz,
                                                      aparc_aseg_in_d_nii_gz[-1],
                                                      "aparc%s_aseg_in_t1_in_b0" % atlas_suffix)

        return job3, jobs4

    def _add_fs_steps(self, dax, job_b0, job_t1, jobs_aparc_aseg):

        b0_nii_gz = File(DWIFiles.B0_NII_GZ.value)
        b0_in_t1_mgz = File(CoregFiles.B0_IN_T1_MGZ.value)
        d2t_reg = File("d2t.reg")
        d2t_lta = File("d2t.lta")
        d2t_mat = File("d2t.mat")
        job1 = Job(CoregJobNames.BBREGISTER.value)
        job1.addArguments(self.subject, b0_nii_gz, b0_in_t1_mgz, d2t_reg, d2t_lta, d2t_mat)
        job1.uses(b0_nii_gz, link=Link.INPUT)
        job1.uses(b0_in_t1_mgz, link=Link.OUTPUT, transfer=False, register=False)
        job1.uses(d2t_reg, link=Link.OUTPUT, transfer=False, register=False)
        job1.uses(d2t_lta, link=Link.OUTPUT, transfer=False, register=False)
        job1.uses(d2t_mat, link=Link.OUTPUT, transfer=False, register=False)
        dax.addJob(job1)

        dax.depends(job1, job_b0)

        b0_in_t1_nii_gz = File(CoregFiles.B0_IN_T1.value)
        job2 = Job(T1JobNames.MRI_CONVERT.value)
        job2.addArguments(b0_in_t1_mgz, b0_in_t1_nii_gz, "--out_orientation", "RAS")
        job2.uses(b0_in_t1_mgz, link=Link.INPUT)
        job2.uses(b0_in_t1_nii_gz, link=Link.OUTPUT, transfer=False, register=False)
        dax.addJob(job2)

        dax.depends(job2, job1)

        # self.qc_snapshots.add_2vols_snapshot_step(dax, [job_t1, job1], t1_nii_gz, b0_in_t1_nii_gz, "b0_in_t1")

        t1_mgz = File(T1Files.T1_MGZ.value)
        t1_in_d_nii_gz = File(CoregFiles.T1_IN_D.value)
        t1_in_d_lta = File(CoregFiles.T1_IN_D.value + ".lta")
        job3 = Job(CoregJobNames.MRI_VOL2VOL.value)
        job3.addArguments("--mov", t1_mgz, "--targ", b0_nii_gz, "--o", t1_in_d_nii_gz, "--lta-inv", d2t_lta,
                          "--save-reg")
        job3.uses(t1_mgz, link=Link.INPUT)
        job3.uses(b0_nii_gz, link=Link.INPUT)
        job3.uses(d2t_lta, link=Link.INPUT)
        job3.uses(t1_in_d_lta, link=Link.OUTPUT, transfer=False, register=False)
        job3.uses(t1_in_d_nii_gz, link=Link.OUTPUT, transfer=False, register=False)
        dax.addJob(job3)

        dax.depends(job3)
        dax.depends(job3, job2)

        # self.qc_snapshots.add_2vols_snapshot_step(dax, [job3], b0_nii_gz, t1_in_d_nii_gz, "t1_in_b0")

        jobs4 = []
        aparc_aseg_mgz = []
        aparc_aseg_in_d_nii_gz = []
        for atlas_suffix, job_aparc_aseg in zip(self.atlas_suffixes, jobs_aparc_aseg):
            aparc_aseg_mgz.append(File(T1Files.APARC_ASEG_MGZ.value % atlas_suffix))
            aparc_aseg_in_d_nii_gz.append(File(CoregFiles.APARC_ASEG_IN_D.value % atlas_suffix))
            jobs4.append(Job(CoregJobNames.MRI_VOL2VOL.value, node_label="Register APARC%s+ASEG to DWI" % atlas_suffix))
            jobs4[-1].addArguments("--mov", aparc_aseg_mgz[-1], "--targ", b0_nii_gz, "--o", aparc_aseg_in_d_nii_gz[-1], "--reg",
                              t1_in_d_lta, "--nearest")
            jobs4[-1].uses(aparc_aseg_mgz[-1], link=Link.INPUT)
            jobs4[-1].uses(b0_nii_gz, link=Link.INPUT)
            jobs4[-1].uses(t1_in_d_lta, link=Link.INPUT)
            jobs4[-1].uses(aparc_aseg_in_d_nii_gz[-1], link=Link.OUTPUT, transfer=False, register=False)
            dax.addJob(jobs4[-1])

            dax.depends(jobs4[-1], job_aparc_aseg)
            dax.depends(jobs4[-1], job3)

            # self.qc_snapshots.add_2vols_snapshot_step(dax, [jobs4[-1]], b0_nii_gz, aparc_aseg_in_d_nii_gz[-1],
            #                                           "aparc%s_aseg_in_b0" % atlas_suffix)
            # self.qc_snapshots.add_3vols_snapshot_step(dax, [job3, jobs4[-1]], t1_in_d_nii_gz, b0_nii_gz,
            #                                           aparc_aseg_in_d_nii_gz[-1], "aparc%s_aseg_in_t1_in_b0" % atlas_suffix)

        return job3, jobs4

    # job_b0 = job6, job_t1 = job7, job_aparc_aseg = job9
    def add_coregistration_steps(self, dax, job_b0, job_t1, job_aparc_aseg):

        if self.use_flirt == "True":
            return self._add_flirt_steps(dax, job_b0, job_t1, job_aparc_aseg)
        else:
            return self._add_fs_steps(dax, job_b0, job_t1, job_aparc_aseg)
