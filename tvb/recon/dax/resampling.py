from Pegasus.DAX3 import File, Link, Job

from tvb.recon.dax.mappings import ResamplingFiles, T1Files, AsegFiles, ResamplingJobNames


class Resampling(object):
    def __init__(self, subject, trg_subject, decim_factor):
        self.subject = subject
        self.trg_subject = trg_subject
        self.decim_factor = decim_factor

    def add_surface_resampling_steps(self, dax, job_recon, job_aseg):
        t1_mgz = File(T1Files.T1_MGZ.value)

        lh_pial = File(T1Files.LH_PIAL.value)
        rh_pial = File(T1Files.RH_PIAL.value)
        pials = [lh_pial, rh_pial]

        # lh_aseg = File(AsegFiles.LH_ASEG.value)
        # rh_aseg = File(AsegFiles.RH_ASEG.value)
        # asegs = [lh_aseg, rh_aseg]

        lh_pial_resamp = File(ResamplingFiles.LH_PIAL_RESAMP.value % self.trg_subject)
        rh_pial_resamp = File(ResamplingFiles.RH_PIAL_RESAMP.value % self.trg_subject)
        pials_resamp = [lh_pial_resamp, rh_pial_resamp]

        # lh_aseg_resamp = File(ResamplingFiles.LH_ASEG_RESAMP.value % self.trg_subject)
        # rh_aseg_resamp = File(ResamplingFiles.RH_ASEG_RESAMP.value % self.trg_subject)
        # asegs_resamp = [lh_aseg_resamp, rh_aseg_resamp]

        for idx, hemi in enumerate(["lh", "rh"]):
            job1 = Job(ResamplingJobNames.MRI_SURF2SURF.value)
            job1.addArguments("--srcsubject", self.subject, "--trgsubject", self.trg_subject, "--hemi", hemi,
                              "--sval-xyz", "pial", "--tval", "pial-%s" % self.trg_subject, "--tval-xyz", t1_mgz)
            job1.uses(t1_mgz, link=Link.INPUT)
            job1.uses(pials[idx], link=Link.INPUT)
            job1.uses(pials_resamp[idx], link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(job1)

            # dax.depends(job1, job_recon)

            # job2 = Job(ResamplingJobNames.MRIS_DECIMATE.value)
            # job2.addArguments("-d", self.decim_factor, asegs[idx], asegs_resamp[idx])
            # job2.uses(asegs[idx], link=Link.INPUT)
            # job2.uses(asegs_resamp[idx], link=Link.OUTPUT, transfer=True, register=True)
            # dax.addJob(job2)
            #
            # dax.depends(job2, job_aseg)

    def add_annotation_resampling_steps(self, dax, job_recon, job_aseg):
        lh_aparc_annot = File(T1Files.LH_APARC_ANNOT.value)
        rh_aparc_annot = File(T1Files.RH_APARC_ANNOT.value)
        aparcs = [lh_aparc_annot, rh_aparc_annot]

        # lh_aseg_annot = File(AsegFiles.LH_ASEG_ANNOT.value)
        # rh_aseg_annot = File(AsegFiles.RH_ASEG_ANNOT.value)
        # asegs = [lh_aseg_annot, rh_aseg_annot]

        lh_aparc_annot_resamp = File(ResamplingFiles.LH_APARC_ANNOT_RESAMP.value % self.trg_subject)
        rh_aparc_annot_resamp = File(ResamplingFiles.RH_APARC_ANNOT_RESAMP.value % self.trg_subject)
        aparcs_resamp = [lh_aparc_annot_resamp, rh_aparc_annot_resamp]

        # lh_aseg_annot_resamp = File(ResamplingFiles.LH_ASEG_ANNOT_RESAMP.value % self.trg_subject)
        # rh_aseg_annot_resamp = File(ResamplingFiles.RH_ASEG_ANNOT_RESAMP.value % self.trg_subject)
        # asegs_resamp = [lh_aseg_annot_resamp, rh_aseg_annot_resamp]

        for idx, hemi in enumerate(["lh", "rh"]):
            job1 = Job(ResamplingJobNames.MRI_SURF2SURF.value)
            job1.addArguments("--srcsubject", self.subject, "--trgsubject", self.trg_subject, "--hemi", hemi,
                              "--sval-annot", aparcs[idx], "--tval", aparcs_resamp[idx])
            job1.uses(aparcs[idx], link=Link.INPUT)
            job1.uses(aparcs_resamp[idx], link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(job1)

            # dax.depends(job1, job_recon)

            # job2 = Job(ResamplingJobNames.?.value)
            # job2.addArguments()
            # job2.uses(asegs[idx], link=Link.INPUT)
            # job2.uses(asegs_resamp[idx], link=Link.OUTPUT, transfer=True, register=True)
            # dax.addJob(job2)
            #
            # dax.depends(job2, job_aseg)