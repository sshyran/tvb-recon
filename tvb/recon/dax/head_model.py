from Pegasus.DAX3 import Job, File, Link

from tvb.recon.dax.mappings import HeadModelJobNames, HeadModelFiles


class HeadModel(object):
    def __init__(self, subject):
        self.subject = subject

    def add_head_model_steps(self, dax, job_recon):
        brain_surface = File(HeadModelFiles.BRAIN_SURFACE.value % self.subject)
        inner_skull_surface = File(HeadModelFiles.INNER_SKULL_SURFACE.value % self.subject)
        outer_skin_surface = File(HeadModelFiles.OUTER_SKIN_SURFACE.value % self.subject)
        outer_skull_surface = File(HeadModelFiles.OUTER_SKULL_SURFACE.value % self.subject)

        bem_surfs = [brain_surface, inner_skull_surface, outer_skin_surface, outer_skull_surface]

        job1 = Job(HeadModelJobNames.MNE_WATERSHED_BEM.value)
        job1.addArguments(self.subject)
        for surf in bem_surfs:
            job1.uses(surf, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job1)

        dax.depends(job1, job_recon)

        brain_surface_tri = File(HeadModelFiles.BRAIN_SURFACE_TRI.value % self.subject)
        inner_skull_surface_tri = File(HeadModelFiles.INNER_SKULL_SURFACE_TRI.value % self.subject)
        outer_skin_surface_tri = File(HeadModelFiles.OUTER_SKIN_SURFACE_TRI.value % self.subject)
        outer_skull_surface_tri = File(HeadModelFiles.OUTER_SKULL_SURFACE_TRI.value % self.subject)

        bem_tri_surfs = [brain_surface_tri, inner_skull_surface_tri, outer_skull_surface_tri, outer_skin_surface_tri]

        last_job = job1

        for idx, surf in enumerate(bem_surfs):
            job = Job(HeadModelJobNames.CONVERT_TO_BRAIN_VISA.value)
            job.addArguments(surf, bem_tri_surfs[idx])
            job.uses(surf, link=Link.INPUT)
            job.uses(bem_tri_surfs[idx], link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(job)

            dax.depends(job, job1)
            last_job = job
