import os
import sys
import time
from Pegasus.DAX3 import ADAG
from tvb.recon.dax import AtlasSuffix, Atlas
from tvb.recon.dax.aseg_generation import AsegGeneration
from tvb.recon.dax.configuration import Configuration, ConfigKey, SensorsType
from tvb.recon.dax.coregistration import Coregistration
from tvb.recon.dax.dwi_processing import DWIProcessing
from tvb.recon.dax.head_model import HeadModel
from tvb.recon.dax.lead_field_model import LeadFieldModel
from tvb.recon.dax.output_conversion import OutputConversion
from tvb.recon.dax.projection_computation import ProjectionComputation
from tvb.recon.dax.resampling import Resampling
from tvb.recon.dax.seeg_computation import SEEGComputation
from tvb.recon.dax.seeg_gain_computation import SeegGainComputation
from tvb.recon.dax.sensor_model import SensorModel
from tvb.recon.dax.source_model import SourceModel
from tvb.recon.dax.t1_processing import T1Processing
from tvb.recon.dax.tracts_generation import TractsGeneration

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: %s DAXFILE\n" % (sys.argv[0]))
        sys.exit(1)
    daxfile = sys.argv[1]
    patient_file = sys.argv[2]

    dax = ADAG("TVB-PIPELINE")
    dax.metadata("created", time.ctime())

    config = Configuration(patient_file)

    subject = config.props[ConfigKey.SUBJECT]
    trg_subject = config.props[ConfigKey.TRGSUBJECT]

    atlas_suffix = AtlasSuffix.DEFAULT

    if config.props[ConfigKey.ATLAS] == Atlas.A2009S:
        atlas_suffix = AtlasSuffix.A2009S

    t1_processing = T1Processing(subject, config.props[ConfigKey.T1_FRMT], config.props[ConfigKey.T2_FLAG],
                                 config.props[ConfigKey.T2_FRMT], config.props[ConfigKey.FLAIR_FLAG],
                                 config.props[ConfigKey.FLAIR_FRMT], config.props[ConfigKey.OPENMP_THRDS],
                                 atlas_suffix)

    resampling = Resampling(subject, trg_subject, config.props[ConfigKey.DECIM_FACTOR], atlas_suffix)

    dwi_processing = DWIProcessing(config.props[ConfigKey.DWI_IS_REVERSED], config.props[ConfigKey.DWI_FRMT],
                                   config.props[ConfigKey.DWI_USE_GRADIENT], config.props[ConfigKey.MRTRIX_THRDS],
                                   config.props[ConfigKey.DWI_SCAN_DIRECTION], config.props[ConfigKey.OS])

    coregistration = Coregistration(subject, config.props[ConfigKey.USE_FLIRT], atlas_suffix)

    tracts_generation = TractsGeneration(config.props[ConfigKey.DWI_MULTI_SHELL], config.props[ConfigKey.MRTRIX_THRDS],
                                         config.props[ConfigKey.STRMLNS_NO], config.props[ConfigKey.STRMLNS_SIFT_NO],
                                         config.props[ConfigKey.STRMLNS_LEN], config.props[ConfigKey.STRMLNS_STEP],
                                         atlas_suffix, config.props[ConfigKey.OS])

    aseg_generation = AsegGeneration(subject, config.props[ConfigKey.ASEG_LH_LABELS],
                                     config.props[ConfigKey.ASEG_RH_LABELS], trg_subject, atlas_suffix)

    output_conversion = OutputConversion(atlas_suffix)

    seeg_computation = SEEGComputation(subject, config.props[ConfigKey.CT_FRMT],
                                       config.props[ConfigKey.CT_ELEC_INTENSITY_TH])

    job_t1, job_aparc_aseg = t1_processing.add_t1_processing_steps(dax, config.props[ConfigKey.RESAMPLE_FLAG])
    job_b0, job_mask = dwi_processing.add_dwi_processing_steps(dax)
    job_t1_in_d, job_aparc_aseg_in_d = coregistration.add_coregistration_steps(dax, job_b0, job_t1, job_aparc_aseg)
    job_aseg_lh, job_aseg_rh = aseg_generation.add_aseg_generation_steps(dax, job_aparc_aseg)

    job_resamp_cort = None
    if config.props[ConfigKey.RESAMPLE_FLAG] == "True":
        job_resamp_cort = resampling.add_surface_resampling_steps(dax, job_aparc_aseg)

    job_mapping_details = aseg_generation.add_mapping_details_computation_step(dax, job_aseg_lh, job_aseg_rh,
                                                                               job_resamp_cort)

    job_weights, job_lengths = tracts_generation.add_tracts_generation_steps(dax, job_t1_in_d, job_mask,
                                                                             job_aparc_aseg_in_d, job_mapping_details)
    job_conn = output_conversion.add_conversion_steps(dax, job_aparc_aseg, job_mapping_details, job_weights,
                                                      job_lengths)

    head_model = None
    job_bem_surfaces = None
    if config.props[ConfigKey.BEM_SURFACES] == "True" or config.props[ConfigKey.USE_OPENMEEG] == "True":
        head_model = HeadModel(subject)
        job_bem_surfaces = head_model.generate_bem_surfaces(dax, job_t1)

    if config.props[ConfigKey.CT_FLAG] == "True":
        if config.props[ConfigKey.SEEG_FLAG] == "True":
            job_seeg_xyz = seeg_computation.add_seeg_positions_computation_steps(dax)

            if config.props[ConfigKey.USE_OPENMEEG] == "True":
                job_head_model = head_model.add_head_model_steps(dax, job_bem_surfaces)

                source_model = SourceModel(subject, trg_subject, atlas_suffix)
                job_source_model = source_model.add_source_model_steps(dax, job_head_model, job_mapping_details)

                sensor_model = SensorModel(subject, trg_subject, atlas_suffix)
                job_sensor_model_lh, job_sensor_model_rh = sensor_model.add_sensor_model_steps(dax, job_source_model)

                lead_field_model = LeadFieldModel(subject, trg_subject, atlas_suffix)
                lead_field_model.add_lead_field_model_steps(dax, job_sensor_model_lh, job_sensor_model_rh)

            else:
                seeg_gain_computation = SeegGainComputation(config.props[ConfigKey.SUBJECT], atlas_suffix)
                if config.props[ConfigKey.SEEG_GAIN_USE_DP] == "True":
                    seeg_gain_computation.add_seeg_gain_dp_computation_steps(dax, job_seeg_xyz, job_mapping_details)
                if config.props[ConfigKey.SEEG_GAIN_USE_MRS] == "True":
                    seeg_gain_computation.add_seeg_mrs_gain_computation_steps(dax, job_seeg_xyz, job_mapping_details)
        else:
            if config.props[ConfigKey.EEG_FLAG] == "True":
                projection_computation = ProjectionComputation(config.props[ConfigKey.SUBJECT], SensorsType.EEG.value,
                                                               atlas_suffix)
                projection_computation.add_projection_computation_steps(dax, job_mapping_details)

            if config.props[ConfigKey.MEG_FLAG] == "True":
                projection_computation = ProjectionComputation(config.props[ConfigKey.SUBJECT], SensorsType.MEG.value,
                                                               atlas_suffix)
                projection_computation.add_projection_computation_steps(dax, job_mapping_details)

    out_dir = os.path.dirname(daxfile)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    with open(daxfile, "w") as f:
        dax.writeXML(f)
