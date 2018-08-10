from Pegasus.DAX3 import Job, File, Link
from tvb.recon.dax.mappings import LeadFieldModelFiles, LeadFieldModelJobNames, HeadModelFiles, SensorModelFiles, \
    SourceModelFiles


class LeadFieldModel(object):
    def __init__(self, subject, trg_subject, atlas_suffix):
        self.subject = subject
        self.trg_subject = trg_subject
        self.atlas_suffix = atlas_suffix
        if len(atlas_suffix) == 0:
            self.atlas_name = "default"
        else:
            self.atlas_name = atlas_suffix[1:]

    def add_lead_field_model_steps(self, dax, job_sensor_model_lh, job_sensor_model_rh):


        head_inv_matrix = File(HeadModelFiles.HEAD_INV_MAT.value)
        head2ipm_file = File(SensorModelFiles.SEEG_H2IPM.value)

        lh_white_dsm = File(SourceModelFiles.LH_WHITE_RESAMP_DSM.value % (self.trg_subject, self.atlas_suffix))
        lh_ds2ipm_file = File(SensorModelFiles.LH_DS2IPM.value % (self.trg_subject, self.atlas_suffix))

        lh_cortical_gain = File(LeadFieldModelFiles.LH_CORT_GAIN_H5.value % self.atlas_suffix)

        job1 = Job(LeadFieldModelJobNames.OM_GAIN.value,
                   node_label="om_gain lh cortical seeg lead field for atlas %s" % self.atlas_name)
        job1.addArguments("-InternalPotential", head_inv_matrix, lh_white_dsm, head2ipm_file, lh_ds2ipm_file,
                          lh_cortical_gain)
        job1.uses(head_inv_matrix, link=Link.INPUT)
        job1.uses(lh_white_dsm, link=Link.INPUT)
        job1.uses(head2ipm_file, link=Link.INPUT)
        job1.uses(lh_ds2ipm_file, link=Link.INPUT)
        job1.uses(lh_cortical_gain, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job1)

        dax.depends(job1, job_sensor_model_lh)

        rh_white_dsm = File(SourceModelFiles.RH_WHITE_RESAMP_DSM.value % (self.trg_subject, self.atlas_suffix))
        rh_ds2ipm_file = File(SensorModelFiles.RH_DS2IPM.value % (self.trg_subject, self.atlas_suffix))

        rh_cortical_gain = File(LeadFieldModelFiles.RH_CORT_GAIN_H5.value % self.atlas_suffix)

        job2 = Job(LeadFieldModelJobNames.OM_GAIN.value,
                   node_label="om_gain rh cortical seeg lead field for atlas %s" % self.atlas_name)
        job2.addArguments("-InternalPotential", head_inv_matrix, rh_white_dsm, head2ipm_file, rh_ds2ipm_file,
                          rh_cortical_gain)
        job2.uses(head_inv_matrix, link=Link.INPUT)
        job2.uses(rh_white_dsm, link=Link.INPUT)
        job2.uses(head2ipm_file, link=Link.INPUT)
        job2.uses(rh_ds2ipm_file, link=Link.INPUT)
        job2.uses(rh_cortical_gain, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job2)

        dax.depends(job2, job_sensor_model_rh)
