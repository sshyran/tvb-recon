from Pegasus.DAX3 import Job, File, Link
from tvb.recon.dax.mappings import SEEGCompFiles, SensorModelFiles, HeadModelJobNames, HeadModelFiles, SourceModelFiles


class SensorModel(object):
    def __init__(self, subject, trg_subject, atlas_suffix):
        self.subject = subject
        self.trg_subject = trg_subject
        self.atlas_suffix = atlas_suffix

    def add_sensor_model_steps(self, dax, job_source_model):
        # TODO: seeg positions file should contain only positions, not labels in order to work with OpenMEEG
        seeg_xyz = File(SEEGCompFiles.SEEG_XYZ.value)

        head_model_geom = File(HeadModelFiles.HEAD_MODEL_GEOM.value)
        head_model_cond = File(HeadModelFiles.HEAD_MODEL_COND.value)

        bem_tri_surfs = [File(HeadModelFiles.INNER_SKULL_SURFACE_LOW_TRI.value % self.subject),
                         File(HeadModelFiles.OUTER_SKULL_SURFACE_LOW_TRI.value % self.subject),
                         File(HeadModelFiles.OUTER_SKIN_SURFACE_LOW_TRI.value % self.subject),
                         File(HeadModelFiles.BRAIN_SURFACE_LOW_TRI.value % self.subject)]

        head2ipm_file = File(SensorModelFiles.SEEG_H2IPM.value)

        job1 = Job(HeadModelJobNames.OM_ASSEMBLE.value)
        job1.addArguments("-h2ipm", head_model_geom, head_model_cond, seeg_xyz, head2ipm_file)
        for surf in bem_tri_surfs:
            job1.uses(surf, link=Link.INPUT)
        job1.uses(head_model_geom, link=Link.INPUT)
        job1.uses(head_model_cond, link=Link.INPUT)
        job1.uses(seeg_xyz, link=Link.INPUT)
        job1.uses(head2ipm_file, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job1)

        dax.depends(job1, job_source_model)

        lh_white_dsm = File(SourceModelFiles.LH_WHITE_RESAMP_DSM.value % (self.trg_subject, self.atlas_suffix))
        lh_ds2ipm_file = File(SensorModelFiles.LH_DS2IPM.value % (self.trg_subject, self.atlas_suffix))

        job2 = Job(HeadModelJobNames.OM_ASSEMBLE.value)
        job2.addArguments("-ds2ipm", head_model_geom, head_model_cond, lh_white_dsm, seeg_xyz, lh_ds2ipm_file)
        for surf in bem_tri_surfs:
            job2.uses(surf, link=Link.INPUT)
        job2.uses(head_model_geom, link=Link.INPUT)
        job2.uses(head_model_cond, link=Link.INPUT)
        job2.uses(lh_white_dsm, link=Link.INPUT)
        job2.uses(seeg_xyz, link=Link.INPUT)
        job2.uses(lh_ds2ipm_file, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job2)

        dax.depends(job2, job1)

        rh_white_dsm = File(SourceModelFiles.RH_WHITE_RESAMP_DSM.value % (self.trg_subject, self.atlas_suffix))
        rh_ds2ipm_file = File(SensorModelFiles.RH_DS2IPM.value % (self.trg_subject, self.atlas_suffix))

        job3 = Job(HeadModelJobNames.OM_ASSEMBLE.value)
        job3.addArguments("-ds2ipm", head_model_geom, head_model_cond, rh_white_dsm, seeg_xyz, rh_ds2ipm_file)
        for surf in bem_tri_surfs:
            job3.uses(surf, link=Link.INPUT)
        job3.uses(head_model_geom, link=Link.INPUT)
        job3.uses(head_model_cond, link=Link.INPUT)
        job3.uses(rh_white_dsm, link=Link.INPUT)
        job3.uses(seeg_xyz, link=Link.INPUT)
        job3.uses(rh_ds2ipm_file, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job3)

        dax.depends(job3, job1)

        return job2, job3
