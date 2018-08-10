from Pegasus.DAX3 import Job, File, Link
from tvb.recon.dax import AtlasSuffix
from tvb.recon.dax.mappings import SEEGCompFiles, SensorModelFiles, HeadModelJobNames, HeadModelFiles, SourceModelFiles


class SensorModel(object):
    def __init__(self, subject, trg_subject, atlas_suffixes=[AtlasSuffix.DEFAULT]):
        self.subject = subject
        self.trg_subject = trg_subject
        self.atlas_suffixes = atlas_suffixes

    def add_sensor_model_steps(self, dax, job_head_model, jobs_source_model):
        # TODO: seeg positions file should contain only positions, not labels in order to work with OpenMEEG
        seeg_xyz = File(SEEGCompFiles.SEEG_XYZ.value)

        head_model_geom = File(HeadModelFiles.HEAD_MODEL_GEOM.value)
        head_model_cond = File(HeadModelFiles.HEAD_MODEL_COND.value)

        bem_tri_surfs = [File(HeadModelFiles.INNER_SKULL_SURFACE_LOW_TRI.value % self.subject),
                         File(HeadModelFiles.OUTER_SKULL_SURFACE_LOW_TRI.value % self.subject),
                         File(HeadModelFiles.OUTER_SKIN_SURFACE_LOW_TRI.value % self.subject),
                         File(HeadModelFiles.BRAIN_SURFACE_LOW_TRI.value % self.subject)]

        head2ipm_file = File(SensorModelFiles.SEEG_H2IPM.value)

        job1 = Job(HeadModelJobNames.OM_ASSEMBLE.value, node_label="om_assemble subcortical seeg sensor model")
        job1.addArguments("-h2ipm", head_model_geom, head_model_cond, seeg_xyz, head2ipm_file)
        for surf in bem_tri_surfs:
            job1.uses(surf, link=Link.INPUT)
        job1.uses(head_model_geom, link=Link.INPUT)
        job1.uses(head_model_cond, link=Link.INPUT)
        job1.uses(seeg_xyz, link=Link.INPUT)
        job1.uses(head2ipm_file, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job1)

        dax.depends(job1, job_head_model)

        lh_white_dsms = []
        rh_white_dsms = []
        lh_ds2ipm_files = []
        rh_ds2ipm_files = []
        jobs_lh = []
        jobs_rh = []
        for iatlas, atlas_suffix in self.atlas_suffixes:

            if len(atlas_suffix) == 0:
                atlas_name = "default"
            else:
                atlas_name = atlas_suffix[1:]

            lh_white_dsms.append(File(SourceModelFiles.LH_WHITE_RESAMP_DSM.value % (self.trg_subject, atlas_suffix)))
            lh_ds2ipm_files.append(File(SensorModelFiles.LH_DS2IPM.value % (self.trg_subject, atlas_suffix)))
    
            jobs_lh.append(Job(HeadModelJobNames.OM_ASSEMBLE.value,
                           node_label="om_assemble lh cortical seeg sensor model for atlas %s" % atlas_name))
            jobs_lh[-1].addArguments("-ds2ipm", head_model_geom, head_model_cond, 
                                     lh_white_dsms[-1], seeg_xyz, lh_ds2ipm_files[-1])
            for surf in bem_tri_surfs:
                jobs_lh[-1].uses(surf, link=Link.INPUT)
            jobs_lh[-1].uses(head_model_geom, link=Link.INPUT)
            jobs_lh[-1].uses(head_model_cond, link=Link.INPUT)
            jobs_lh[-1].uses(lh_white_dsms[-1], link=Link.INPUT)
            jobs_lh[-1].uses(seeg_xyz, link=Link.INPUT)
            jobs_lh[-1].uses(lh_ds2ipm_files[-1], link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(jobs_lh[-1])
    
            dax.depends(jobs_lh[-1], job1)
            dax.depends(jobs_lh[-1], jobs_source_model[0][iatlas])
    
            rh_white_dsms.append(File(SourceModelFiles.RH_WHITE_RESAMP_DSM.value % (self.trg_subject, atlas_suffix)))
            rh_ds2ipm_files.append(File(SensorModelFiles.RH_DS2IPM.value % (self.trg_subject, atlas_suffix)))
    
            jobs_rh.append(Job(HeadModelJobNames.OM_ASSEMBLE.value,
                           node_label="om_assemble rh cortical seeg sensor model for atlas %s" % atlas_name))
            jobs_rh[-1].addArguments("-ds2ipm", head_model_geom, head_model_cond, rh_white_dsms[-1], seeg_xyz, rh_ds2ipm_files[-1])
            for surf in bem_tri_surfs:
                jobs_rh[-1].uses(surf, link=Link.INPUT)
            jobs_rh[-1].uses(head_model_geom, link=Link.INPUT)
            jobs_rh[-1].uses(head_model_cond, link=Link.INPUT)
            jobs_rh[-1].uses(rh_white_dsms[-1], link=Link.INPUT)
            jobs_rh[-1].uses(seeg_xyz, link=Link.INPUT)
            jobs_rh[-1].uses(rh_ds2ipm_files[-1], link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(jobs_rh[-1])
    
            dax.depends(jobs_rh[-1], job1)
            dax.depends(jobs_lh[-1], jobs_source_model[1][iatlas])

        return jobs_lh, jobs_rh
