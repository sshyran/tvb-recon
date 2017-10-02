import os
import sys
import time
from Pegasus.DAX3 import ADAG
from tvb.recon.dax.aseg_generation import AsegGeneration
from tvb.recon.dax.configuration import Configuration, ConfigKey, SensorsType
from tvb.recon.dax.coregistration import Coregistration
from tvb.recon.dax.dwi_processing import DWIProcessing
from tvb.recon.dax.head_model import HeadModel
from tvb.recon.dax.output_conversion import OutputConversion
from tvb.recon.dax.projection_computation import ProjectionComputation
from tvb.recon.dax.resampling import Resampling
from tvb.recon.dax.seeg_computation import SEEGComputation
from tvb.recon.dax.seeg_gain_computation import SeegGainComputation
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

    t1_processing = T1Processing(config.props[ConfigKey.SUBJECT], config.props[ConfigKey.T1_FRMT],
                                 config.props[ConfigKey.T2_FLAG], config.props[ConfigKey.T2_FRMT],
                                 config.props[ConfigKey.FLAIR_FLAG], config.props[ConfigKey.FLAIR_FRMT],
                                 config.props[ConfigKey.OPENMP_THRDS])

    resampling = Resampling(config.props[ConfigKey.SUBJECT], config.props[ConfigKey.TRGSUBJECT],
                            config.props[ConfigKey.DECIM_FACTOR])

    dwi_processing = DWIProcessing(config.props[ConfigKey.DWI_IS_REVERSED], config.props[ConfigKey.DWI_FRMT],
                                   config.props[ConfigKey.MRTRIX_THRDS], config.props[ConfigKey.DWI_SCAN_DIRECTION])

    coregistration = Coregistration(config.props[ConfigKey.SUBJECT], config.props[ConfigKey.USE_FLIRT])

    tracts_generation = TractsGeneration(config.props[ConfigKey.DWI_MULTI_SHELL], config.props[ConfigKey.MRTRIX_THRDS],
                                         config.props[ConfigKey.STRMLNS_NO], config.props[ConfigKey.STRMLNS_SIFT_NO],
                                         config.props[ConfigKey.STRMLNS_LEN], config.props[ConfigKey.STRMLNS_STEP])

    aseg_generation = AsegGeneration(config.props[ConfigKey.SUBJECT], config.props[ConfigKey.ASEG_LH_LABELS],
                                     config.props[ConfigKey.ASEG_RH_LABELS])

    output_conversion = OutputConversion()

    seeg_computation = SEEGComputation(config.props[ConfigKey.SUBJECT], config.props[ConfigKey.CT_FRMT],
                                       config.props[ConfigKey.CT_ELEC_INTENSITY_TH])

    job_t1, job_aparc_aseg = t1_processing.add_t1_processing_steps(dax)
    job_b0, job_mask = dwi_processing.add_dwi_processing_steps(dax)
    job_t1_in_d, job_aparc_aseg_in_d = coregistration.add_coregistration_steps(dax, job_b0, job_t1,
                                                                               job_aparc_aseg)
    job_aseg_lh, job_aseg_rh, job_mapping_details = aseg_generation.add_aseg_generation_steps(dax, job_aparc_aseg)

    if config.props[ConfigKey.RESAMPLE_FLAG] == "True":
        resampling.add_surface_resampling_steps(dax, job_aparc_aseg, job_aseg_lh)
        resampling.add_annotation_resampling_steps(dax, job_aparc_aseg, job_aseg_lh)

    job_weights, job_lengths = tracts_generation.add_tracts_generation_steps(dax, job_t1_in_d, job_mask,
                                                                             job_aparc_aseg_in_d, job_mapping_details)
    job_conn = output_conversion.add_conversion_steps(dax, job_aparc_aseg, job_mapping_details, job_weights,
                                                      job_lengths)

    if config.props[ConfigKey.CT_FLAG] == "True":
        job_seeg_xyz = seeg_computation.add_seeg_positions_computation_steps(dax)

        if config.props[ConfigKey.USE_OPENMEEG] == "True":
            head_model = HeadModel(config.props[ConfigKey.SUBJECT])
            head_model.add_head_model_steps(dax, job_t1)

        else:
            if config.props[ConfigKey.SEEG_FLAG] == "True":
                seeg_gain_computation = SeegGainComputation(config.props[ConfigKey.SUBJECT])
                seeg_gain_computation.add_seeg_gain_computation_steps(dax, job_seeg_xyz, job_mapping_details)
            if config.props[ConfigKey.EEG_FLAG] == "True":
                projection_computation = ProjectionComputation(config.props[ConfigKey.SUBJECT], SensorsType.EEG.value)
                projection_computation.add_projection_computation_steps(dax, job_mapping_details)

            if config.props[ConfigKey.MEG_FLAG] == "True":
                projection_computation = ProjectionComputation(config.props[ConfigKey.SUBJECT], SensorsType.MEG.value)
                projection_computation.add_projection_computation_steps(dax, job_mapping_details)

    out_dir = os.path.dirname(daxfile)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    with open(daxfile, "w") as f:
        dax.writeXML(f)
