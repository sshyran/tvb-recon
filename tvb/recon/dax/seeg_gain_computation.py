from Pegasus.DAX3 import File, Job, Link

from tvb.recon.dax.mappings import SeegGainFiles, SEEGCompFiles, AsegFiles, ProjectionCompJobNames, SeegGainJobNames


class SeegGainComputation(object):
    def __init__(self, subject):
        self.subject = subject

    def add_seeg_gain_dp_computation_steps(self, dax, job_seeg_xyz, job_mapping_details):
        seeg_xyz = File(SEEGCompFiles.SEEG_XYZ.value)
        centers_txt = File(AsegFiles.CENTERS_TXT.value)

        gain_mat = File(SeegGainFiles.SEEG_GAIN_DP_MAT.value)

        job = Job(ProjectionCompJobNames.COMPUTE_PROJ_MAT.value)
        job.addArguments(seeg_xyz, centers_txt, gain_mat, self.subject)
        job.uses(seeg_xyz, link=Link.INPUT)
        job.uses(centers_txt, link=Link.INPUT)
        job.uses(gain_mat, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job)

        dax.depends(job, job_seeg_xyz)
        dax.depends(job, job_mapping_details)

    def add_seeg_mrs_gain_computation_steps(self, dax, job_seeg_xyz, job_mapping_details):
        seeg_xyz = File(SEEGCompFiles.SEEG_XYZ.value)
        cort_surf = File(AsegFiles.SURF_CORT_ZIP.value)
        subcort_surf = File(AsegFiles.SURF_SUBCORT_ZIP.value)
        cort_rm = File(AsegFiles.RM_CORT_TXT.value)
        subcort_rm = File(AsegFiles.RM_SUBCORT_TXT.value)

        gain_mat = File(SeegGainFiles.SEEG_GAIN_MRS_MAT.value)

        job = Job(SeegGainJobNames.COMPUTE_SEEG_GAIN.value)
        job.addArguments(seeg_xyz, cort_surf, subcort_surf, cort_rm, subcort_rm, gain_mat, self.subject)
        job.uses(seeg_xyz, link=Link.INPUT)
        job.uses(cort_surf, link=Link.INPUT)
        job.uses(subcort_surf, link=Link.INPUT)
        job.uses(cort_rm, link=Link.INPUT)
        job.uses(subcort_rm, link=Link.INPUT)
        job.uses(gain_mat, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job)

        dax.depends(job, job_seeg_xyz)
        dax.depends(job, job_mapping_details)
