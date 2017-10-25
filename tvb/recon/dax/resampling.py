from Pegasus.DAX3 import File, Link, Job

from tvb.recon.dax.mappings import ResamplingFiles, T1Files, ResamplingJobNames, T1JobNames
from tvb.recon.dax.qc_snapshots import QCSnapshots


class Resampling(object):
    def __init__(self, subject, trg_subject, decim_factor):
        self.subject = subject
        self.trg_subject = trg_subject
        self.decim_factor = decim_factor
        self.qc_snapshots = QCSnapshots.get_instance()

    def add_surface_resampling_steps(self, dax, job_recon):
        t1_mgz = File(T1Files.T1_MGZ.value)

        lh_pial = File(T1Files.LH_PIAL.value)
        rh_pial = File(T1Files.RH_PIAL.value)
        pials = [lh_pial, rh_pial]

        lh_aparc_annot = File(T1Files.LH_APARC_ANNOT.value)
        rh_aparc_annot = File(T1Files.RH_APARC_ANNOT.value)
        aparcs = [lh_aparc_annot, rh_aparc_annot]

        lh_pial_resamp = File(ResamplingFiles.LH_PIAL_RESAMP.value % self.trg_subject)
        rh_pial_resamp = File(ResamplingFiles.RH_PIAL_RESAMP.value % self.trg_subject)
        pials_resamp = [lh_pial_resamp, rh_pial_resamp]

        lh_aparc_annot_resamp = File(ResamplingFiles.LH_APARC_ANNOT_RESAMP.value % self.trg_subject)
        rh_aparc_annot_resamp = File(ResamplingFiles.RH_APARC_ANNOT_RESAMP.value % self.trg_subject)
        aparcs_resamp = [lh_aparc_annot_resamp, rh_aparc_annot_resamp]

        last_job = None

        for idx, hemi in enumerate(["lh", "rh"]):
            job1 = Job(ResamplingJobNames.MRI_SURF2SURF.value)
            job1.addArguments("--srcsubject", self.subject, "--trgsubject", self.trg_subject, "--hemi", hemi,
                              "--sval-xyz", "pial", "--tval", "pial-%s" % self.trg_subject, "--tval-xyz", t1_mgz)
            job1.uses(t1_mgz, link=Link.INPUT)
            job1.uses(pials[idx], link=Link.INPUT)
            job1.uses(pials_resamp[idx], link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(job1)

            dax.depends(job1, job_recon)

            job2 = Job(ResamplingJobNames.MRI_SURF2SURF.value)
            job2.addArguments("--srcsubject", self.subject, "--trgsubject", self.trg_subject, "--hemi", hemi,
                              "--sval-annot", aparcs[idx], "--tval", aparcs_resamp[idx])
            job2.uses(aparcs[idx], link=Link.INPUT)
            job2.uses(aparcs_resamp[idx], link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(job2)

            dax.depends(job2, job_recon)
            last_job = job2

        lh_centered_pial = File(ResamplingFiles.LH_CENTERED_PIAL_RESAMP.value % self.trg_subject)
        job5 = Job(T1JobNames.MRIS_CONVERT.value)
        job5.addArguments("--to-scanner", lh_pial_resamp, lh_centered_pial)
        job5.uses(lh_pial_resamp, link=Link.INPUT)
        job5.uses(lh_centered_pial, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job5)

        dax.depends(job5, last_job)

        rh_centered_pial = File(ResamplingFiles.RH_CENTERED_PIAL_RESAMP.value % self.trg_subject)
        job6 = Job(T1JobNames.MRIS_CONVERT.value)
        job6.addArguments("--to-scanner", rh_pial_resamp, rh_centered_pial)
        job6.uses(rh_pial_resamp, link=Link.INPUT)
        job6.uses(rh_centered_pial, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job6)

        dax.depends(job6, last_job)

        t1_nii_gz = File(T1Files.T1_NII_GZ.value)
        self.qc_snapshots.add_vol_surf_snapshot_step(dax, [job5, job6], t1_nii_gz, [lh_centered_pial, rh_centered_pial])
        self.qc_snapshots.add_surf_annot_snapshot_step(dax, [job5, job6], lh_centered_pial, lh_aparc_annot_resamp)
        self.qc_snapshots.add_surf_annot_snapshot_step(dax, [job5, job6], rh_centered_pial, rh_aparc_annot_resamp)

        return job6


        # def add_annotation_resampling_steps(self, dax, job_aseg):
        # lh_aseg = File(AsegFiles.LH_ASEG.value)
        # rh_aseg = File(AsegFiles.RH_ASEG.value)
        # asegs = [lh_aseg, rh_aseg]

        # lh_aseg_annot = File(AsegFiles.LH_ASEG_ANNOT.value)
        # rh_aseg_annot = File(AsegFiles.RH_ASEG_ANNOT.value)
        # asegs = [lh_aseg_annot, rh_aseg_annot]

        # lh_aseg_resamp = File(ResamplingFiles.LH_ASEG_RESAMP.value % self.trg_subject)
        # rh_aseg_resamp = File(ResamplingFiles.RH_ASEG_RESAMP.value % self.trg_subject)
        # asegs_resamp = [lh_aseg_resamp, rh_aseg_resamp]

        # lh_aseg_annot_resamp = File(ResamplingFiles.LH_ASEG_ANNOT_RESAMP.value % self.trg_subject)
        # rh_aseg_annot_resamp = File(ResamplingFiles.RH_ASEG_ANNOT_RESAMP.value % self.trg_subject)
        # asegs_resamp = [lh_aseg_annot_resamp, rh_aseg_annot_resamp]

        # for idx, hemi in enumerate(["lh", "rh"]):
        # job2 = Job(ResamplingJobNames.MRIS_DECIMATE.value)
        # job2.addArguments("-d", self.decim_factor, asegs[idx], asegs_resamp[idx])
        # job2.uses(asegs[idx], link=Link.INPUT)
        # job2.uses(asegs_resamp[idx], link=Link.OUTPUT, transfer=True, register=True)
        # dax.addJob(job2)
        #
        # dax.depends(job2, job_aseg)

        # job2 = Job(ResamplingJobNames.?.value)
        # job2.addArguments()
        # job2.uses(asegs[idx], link=Link.INPUT)
        # job2.uses(asegs_resamp[idx], link=Link.OUTPUT, transfer=True, register=True)
        # dax.addJob(job2)
        #
        # dax.depends(job2, job_aseg)
