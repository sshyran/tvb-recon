from Pegasus.DAX3 import Job, File, Link
from tvb.recon.dax.mappings import HeadModelJobNames, HeadModelFiles, ResamplingJobNames


class HeadModel(object):
    def __init__(self, subject, resample):
        self.subject = subject
        self.resample = resample

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

        brain_surface_low = File(HeadModelFiles.BRAIN_SURFACE_LOW.value % self.subject)
        inner_skull_surface_low = File(HeadModelFiles.INNER_SKULL_SURFACE_LOW.value % self.subject)
        outer_skin_surface_low = File(HeadModelFiles.OUTER_SKIN_SURFACE_LOW.value % self.subject)
        outer_skull_surface_low = File(HeadModelFiles.OUTER_SKULL_SURFACE_LOW.value % self.subject)

        bem_surfs_low = [brain_surface_low, inner_skull_surface_low, outer_skin_surface_low, outer_skull_surface_low]

        last_job = job1

        for i, surf in enumerate(bem_surfs):
            job_resamp = Job(ResamplingJobNames.MRIS_DECIMATE.value)
            job_resamp.addArguments("-d", "0.1", surf, bem_surfs_low[i])
            job_resamp.uses(surf, link=Link.INPUT)
            job_resamp.uses(bem_surfs_low[i], link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(job_resamp)

            dax.depends(job_resamp, job1)
            last_job = job_resamp

        brain_surface_tri = File(HeadModelFiles.BRAIN_SURFACE_LOW_TRI.value % self.subject)
        inner_skull_surface_tri = File(HeadModelFiles.INNER_SKULL_SURFACE_LOW_TRI.value % self.subject)
        outer_skin_surface_tri = File(HeadModelFiles.OUTER_SKIN_SURFACE_LOW_TRI.value % self.subject)
        outer_skull_surface_tri = File(HeadModelFiles.OUTER_SKULL_SURFACE_LOW_TRI.value % self.subject)

        bem_tri_surfs = [brain_surface_tri, inner_skull_surface_tri, outer_skull_surface_tri, outer_skin_surface_tri]

        for idx, surf in enumerate(bem_surfs_low):
            tri_file = bem_tri_surfs[idx]
            job2 = Job(HeadModelJobNames.CONVERT_TO_BRAIN_VISA.value)
            job2.addArguments(surf, tri_file, self.subject)
            job2.uses(surf, link=Link.INPUT)
            job2.uses(tri_file, link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(job2)

            dax.depends(job2, job1)
            last_job = job2

        head_model_geom = File(HeadModelFiles.HEAD_MODEL_GEOM.value)
        head_model_cond = File(HeadModelFiles.HEAD_MODEL_COND.value)
        job4 = Job(HeadModelJobNames.GEN_HEAD_MODEL.value)
        job4.addArguments(self.subject, self.resample)
        for surf in bem_tri_surfs:
            job4.uses(surf, link=Link.INPUT)
        job4.uses(head_model_geom, link=Link.OUTPUT, transfer=True, register=True)
        job4.uses(head_model_cond, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job4)

        dax.depends(job4, last_job)

        head_matrix = File(HeadModelFiles.HEAD_MAT.value)
        job5 = Job(HeadModelJobNames.OM_ASSEMBLE.value)
        job5.addArguments("-HM", head_model_geom, head_model_cond, head_matrix)
        for surf in bem_tri_surfs:
            job5.uses(surf, link=Link.INPUT)
        job5.uses(head_model_geom, link=Link.INPUT)
        job5.uses(head_model_cond, link=Link.INPUT)
        job5.uses(head_matrix, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job5)

        dax.depends(job5, last_job)
        dax.depends(job5, job4)

        head_inv_matrix = File(HeadModelFiles.HEAD_INV_MAT.value)
        job6 = Job(HeadModelJobNames.OM_MINVERSER.value)
        job6.addArguments(head_matrix, head_inv_matrix)
        job6.uses(head_matrix, link=Link.INPUT)
        job6.uses(head_inv_matrix, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job6)

        dax.depends(job6, job5)

        return job6
