from Pegasus.DAX3 import File, Job, Link

from tvb.recon.dax.mappings import T1Files, SeegGainFiles, SeegGainJobNames, SEEGCompFiles


class SeegGainComputation(object):
    def __init__(self, subject):
        self.subject = subject

    def add_seeg_gain_computation_steps(self, dax, job_seeg_xyz):

        lh_pial = File(T1Files.LH_PIAL.value)
        rh_pial = File(T1Files.RH_PIAL.value)
        cortical_pial = File(SeegGainFiles.CORTICAL_PIAL.value)
        job1 = Job(SeegGainJobNames.MERGE_SURFACES.value)
        job1.addArguments(lh_pial, rh_pial, cortical_pial, self.subject)
        job1.uses(lh_pial, link=Link.INPUT)
        job1.uses(rh_pial, link=Link.INPUT)
        job1.uses(cortical_pial, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job1)

        dax.depends(job1, job_seeg_xyz)

        lh_aparc_annot = File(T1Files.LH_APARC_ANNOT.value)
        rh_aparc_annot = File(T1Files.RH_APARC_ANNOT.value)
        seeg_xyz = File(SEEGCompFiles.SEEG_XYZ.value)
        gain_mat = File(SeegGainFiles.SEEG_GAIN_MAT.value)
        job2 = Job(SeegGainJobNames.COMPUTE_SEEG_GAIN.value)
        job2.addArguments(cortical_pial, lh_aparc_annot, rh_aparc_annot, seeg_xyz, gain_mat, self.subject)
        job2.uses(cortical_pial, link=Link.INPUT)
        job2.uses(lh_aparc_annot, link=Link.INPUT)
        job2.uses(rh_aparc_annot, link=Link.INPUT)
        job2.uses(seeg_xyz, link=Link.INPUT)
        job2.uses(gain_mat, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job2)

        dax.depends(job2, job1)