from Pegasus.DAX3 import File, Job, Link

from tvb.recon.dax.mappings import T1Files, SeegGainFiles, SeegGainJobNames, SEEGCompFiles, AsegFiles


class SeegGainComputation(object):
    def __init__(self, subject):
        self.subject = subject

    def add_seeg_gain_computation_steps(self, dax, job_seeg_xyz, job_recon, job_aseg_lh, job_aseg_rh):

        lh_pial = File(T1Files.LH_PIAL.value)
        rh_pial = File(T1Files.RH_PIAL.value)
        cortical_pial = File(SeegGainFiles.CORTICAL_PIAL.value)
        job1 = Job(SeegGainJobNames.MERGE_SURFACES.value)
        job1.addArguments(lh_pial, rh_pial, cortical_pial, self.subject)
        job1.uses(lh_pial, link=Link.INPUT)
        job1.uses(rh_pial, link=Link.INPUT)
        job1.uses(cortical_pial, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job1)

        dax.depends(job1, job_recon)

        lh_aparc_annot = File(T1Files.LH_APARC_ANNOT.value)
        rh_aparc_annot = File(T1Files.RH_APARC_ANNOT.value)
        seeg_xyz = File(SEEGCompFiles.SEEG_XYZ.value)
        gain_mat = File(SeegGainFiles.SEEG_GAIN_MAT.value)
        lh_aseg = File(AsegFiles.LH_ASEG.value)
        rh_aseg = File(AsegFiles.RH_ASEG.value)
        subcortical_pial = File(SeegGainFiles.SUBCORTICAL_ASEG.value)

        job = Job(SeegGainJobNames.MERGE_SURFACES.value)
        job.addArguments(lh_aseg, rh_aseg, subcortical_pial, self.subject)
        job.uses(lh_aseg, link=Link.INPUT)
        job.uses(rh_aseg, link=Link.INPUT)
        job.uses(subcortical_pial, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job)

        dax.depends(job, job_aseg_lh)
        dax.depends(job, job_aseg_rh)

        lh_aseg_annot = File(AsegFiles.LH_ASEG_ANNOT.value)
        rh_aseg_annot = File(AsegFiles.RH_ASEG_ANNOT.value)
        job3 = Job(SeegGainJobNames.COMPUTE_SEEG_GAIN.value)
        job3.addArguments(cortical_pial, lh_aparc_annot, rh_aparc_annot, subcortical_pial, lh_aseg_annot, rh_aseg_annot, seeg_xyz, gain_mat, self.subject)
        job3.uses(cortical_pial, link=Link.INPUT)
        job3.uses(lh_aparc_annot, link=Link.INPUT)
        job3.uses(rh_aparc_annot, link=Link.INPUT)
        job3.uses(subcortical_pial, link=Link.INPUT)
        job3.uses(lh_aseg_annot, link=Link.INPUT)
        job3.uses(rh_aseg_annot, link=Link.INPUT)
        job3.uses(seeg_xyz, link=Link.INPUT)
        job3.uses(gain_mat, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job3)

        dax.depends(job3, job_seeg_xyz)
        dax.depends(job3, job)
