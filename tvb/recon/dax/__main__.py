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
from tvb.recon.dax.mrielec_seeg_computation import MriElecSEEGComputation
from tvb.recon.dax.output_conversion import OutputConversion
from tvb.recon.dax.projection_computation import ProjectionComputation
from tvb.recon.dax.resampling import Resampling
from tvb.recon.dax.seeg_computation import SEEGComputation
from tvb.recon.dax.seeg_gain_computation import SeegGainComputation
from tvb.recon.dax.sensor_model import SensorModel
from tvb.recon.dax.source_model import SourceModel
from tvb.recon.dax.t1_processing import T1Processing
from tvb.recon.dax.tracts_generation import TractsGeneration


def get_atlases_and_suffixes_lists(atlases):
    atlases_list = atlases.split(" ")
    atlas_suffixes = []
    for atlas in atlases_list:
        if atlas == Atlas.A2009S:
            atlas_suffixes.append(AtlasSuffix.A2009S)
        elif atlas == Atlas.DKT:
            atlas_suffixes.append(AtlasSuffix.DKT)
        else:
            atlas_suffixes.append(AtlasSuffix.DEFAULT)
    return atlases_list, atlas_suffixes


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

    atlases = ConfigKey.ATLAS
    atlases_list, atlas_suffixes = get_atlases_and_suffixes_lists(atlases)

    #---------------------First we add jobs to the dax for all desired atlases together---------------------------------

    t1_processing = T1Processing(subject, config.props[ConfigKey.T1_FRMT], config.props[ConfigKey.T2_FLAG],
                                 config.props[ConfigKey.T2_FRMT], config.props[ConfigKey.FLAIR_FLAG],
                                 config.props[ConfigKey.FLAIR_FRMT], config.props[ConfigKey.OPENMP_THRDS],
                                 atlas_suffixes)

    job_t1, jobs_aparc_aseg = t1_processing.add_t1_processing_steps(dax, config.props[ConfigKey.RESAMPLE_FLAG])

    dwi_processing = DWIProcessing(config.props[ConfigKey.DWI_IS_REVERSED], config.props[ConfigKey.DWI_FRMT],
                                   config.props[ConfigKey.DWI_USE_GRADIENT], config.props[ConfigKey.MRTRIX_THRDS],
                                   config.props[ConfigKey.DWI_SCAN_DIRECTION], config.props[ConfigKey.OS])
    job_b0, job_mask = dwi_processing.add_dwi_processing_steps(dax)

    coregistration = Coregistration(subject, config.props[ConfigKey.USE_FLIRT], atlas_suffixes)
    job_t1_in_d, jobs_aparc_aseg_in_d = coregistration.add_coregistration_steps(dax, job_b0, job_t1, jobs_aparc_aseg)

    #--------------From this point and on we should add jobs to the dax for each desired atlas separetely---------------

    if config.props[ConfigKey.RESAMPLE_FLAG] == "True":
        job_resamp_cort = []
        resampling = []
    else:
        job_resamp_cort = [None]
        resampling = None

    aseg_generation = []
    job_aseg_lh = []
    job_aseg_rh = []
    job_mapping_details = []
    tracts_generation = []
    job_weights = []
    job_lengths = []
    job_conn = []
    output_conversion = []

    mrielec_seeg_computation = MriElecSEEGComputation(subject, config.props[ConfigKey.MRIELEC_FRMT],
                                                      config.props[ConfigKey.SAME_SPACE_VOL_POM],
                                                      config.props[ConfigKey.CT_ELEC_INTENSITY_TH])

    seeg_computation = SEEGComputation(subject, config.props[ConfigKey.CT_FRMT],
                                       config.props[ConfigKey.CT_ELEC_INTENSITY_TH])

    head_model = None
    job_bem_surfaces = None
    job_source_model = None
    job_sensor_model_lh = None
    job_sensor_model_rh = None
    job_seeg_xyz = None
    projection_computation_EEG = None
    projection_computation_MEG = None

    if config.props[ConfigKey.BEM_SURFACES] == "True" or config.props[ConfigKey.USE_OPENMEEG] == "True":
        head_model = HeadModel(subject)
        job_bem_surfaces = head_model.generate_bem_surfaces(dax, job_t1)
        
    if config.props[ConfigKey.PROCESS_SENSORS] == "True":
        if config.props[ConfigKey.SEEG_FLAG] == "True":
            if config.props[ConfigKey.USE_OPENMEEG] == "True":
                job_head_model = head_model.add_head_model_steps(dax, job_bem_surfaces)
                source_model = []
                job_source_model = []
                sensor_model = []
                job_sensor_model_lh = []
                job_sensor_model_rh = []
                lead_field_model = []
            else:
                seeg_gain_computation = []
                if config.props[ConfigKey.MRIELEC_FLAG] == "True":
                    job_seeg_xyz = mrielec_seeg_computation.add_computation_steps(dax)
                else:
                    if config.props[ConfigKey.CT_FLAG] == "True":
                        job_seeg_xyz = seeg_computation.add_seeg_positions_computation_steps(dax)
        else:
            if config.props[ConfigKey.EEG_FLAG] == "True":
                projection_computation_EEG = []
            if config.props[ConfigKey.MEG_FLAG] == "True":
                projection_computation_MEG = []

    # Having prepared the necessary job lists, loop now for each desired atlas:
    # Always create new objects to make sure we don't create a mess...
    for atlas_suffix, job_aparc_aseg, job_aparc_aseg_in_d in \
        zip(atlas_suffixes, jobs_aparc_aseg, jobs_aparc_aseg_in_d):

        if config.props[ConfigKey.RESAMPLE_FLAG] == "True":
            resampling.append(Resampling(subject, trg_subject, config.props[ConfigKey.DECIM_FACTOR], atlas_suffix))
            job_resamp_cort.append(resampling[-1].add_surface_resampling_steps(dax, job_aparc_aseg))

        aseg_generation.append(AsegGeneration(subject, config.props[ConfigKey.ASEG_LH_LABELS],
                                         config.props[ConfigKey.ASEG_RH_LABELS], trg_subject, atlas_suffix))
        temp_job_aseg_lh, temp_job_aseg_rh = aseg_generation[-1].add_aseg_generation_steps(dax, job_aparc_aseg)
        job_aseg_lh.append(temp_job_aseg_lh)
        job_aseg_rh.append(temp_job_aseg_rh)

        job_mapping_details.append(aseg_generation[-1]. \
                                    add_mapping_details_computation_step(dax, job_aseg_lh[-1], job_aseg_rh[-1],
                                                                         job_resamp_cort[-1]))

        tracts_generation.append(TractsGeneration(config.props[ConfigKey.DWI_MULTI_SHELL],
                                                  config.props[ConfigKey.MRTRIX_THRDS],
                                                  config.props[ConfigKey.STRMLNS_NO],
                                                  config.props[ConfigKey.STRMLNS_SIFT_NO],
                                                  config.props[ConfigKey.STRMLNS_LEN],
                                                  config.props[ConfigKey.STRMLNS_STEP],
                                                  atlas_suffix, config.props[ConfigKey.OS]))

        temp_job_weights, temp_job_lengths = \
            tracts_generation[-1].add_tracts_generation_steps(dax, job_t1_in_d, job_mask,
                                                              job_aparc_aseg_in_d, job_mapping_details[-1])
        job_weights.append(temp_job_weights)
        job_lengths.append(temp_job_lengths)

        output_conversion.append(OutputConversion(atlas_suffix))
        job_conn.append(output_conversion[-1].add_conversion_steps(dax, job_aparc_aseg, job_mapping_details[-1],
                                                                   job_weights[-1], job_lengths[-1]))

        if config.props[ConfigKey.PROCESS_SENSORS] == "True":
            if config.props[ConfigKey.SEEG_FLAG] == "True":
                if config.props[ConfigKey.USE_OPENMEEG] == "True":

                    source_model.append(SourceModel(subject, trg_subject, atlas_suffix))
                    job_source_model.append(source_model[-1].add_source_model_steps(dax, job_head_model,
                                                                                job_mapping_details[-1]))

                    sensor_model.append(SensorModel(subject, trg_subject, atlas_suffix))
                    temp_job_sensor_model_lh, temp_job_sensor_model_rh = \
                        sensor_model[-1].add_sensor_model_steps(dax, job_source_model[-1])
                    job_sensor_model_lh.append(temp_job_sensor_model_lh)
                    job_sensor_model_rh.append(temp_job_sensor_model_rh)

                    lead_field_model.append(LeadFieldModel(subject, trg_subject, atlas_suffix))
                    lead_field_model[-1].add_lead_field_model_steps(dax, job_sensor_model_lh[-1], job_sensor_model_rh[-1])

                else:
                    if job_seeg_xyz is not None:
                        seeg_gain_computation.append(SeegGainComputation(config.props[ConfigKey.SUBJECT], atlas_suffix))
                        if config.props[ConfigKey.SEEG_GAIN_USE_DP] == "True":
                            seeg_gain_computation[-1].add_seeg_gain_dp_computation_steps(dax, job_seeg_xyz,
                                                                                     job_mapping_details[-1])
                        if config.props[ConfigKey.SEEG_GAIN_USE_MRS] == "True":
                            seeg_gain_computation[-1].add_seeg_mrs_gain_computation_steps(dax, job_seeg_xyz,
                                                                                      job_mapping_details[-1])

            else:

                if config.props[ConfigKey.EEG_FLAG] == "True":
                    projection_computation_EEG.append(ProjectionComputation(config.props[ConfigKey.SUBJECT],
                                                                            SensorsType.EEG.value,
                                                                            atlas_suffix))
                    projection_computation_EEG[-1].add_projection_computation_steps(dax, job_mapping_details[-1])

                if config.props[ConfigKey.MEG_FLAG] == "True":
                    projection_computation_MEG.append(ProjectionComputation(config.props[ConfigKey.SUBJECT],
                                                                            SensorsType.MEG.value,
                                                                            atlas_suffix))
                    projection_computation_MEG[-1].add_projection_computation_steps(dax, job_mapping_details[-1])

    out_dir = os.path.dirname(daxfile)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    with open(daxfile, "w") as f:
        dax.writeXML(f)
