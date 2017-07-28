import sys

import time
from Pegasus.DAX3 import ADAG
from t1_processing import T1Processing
from dwi_processing import DWIProcessing
from coregistration import Coregistration
from configuration import Configuration, ConfigKey
from tracts_generation import TractsGeneration

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: %s DAXFILE\n" % (sys.argv[0]))
        sys.exit(1)
    daxfile = sys.argv[1]
    patient_file = sys.argv[2]

    dax = ADAG("convert")
    dax.metadata("created", time.ctime())

    config = Configuration(patient_file)

    t1_processing = T1Processing(t1_frmt=config.props[ConfigKey.T1_FRMT], use_t2=config.props[ConfigKey.T2_FLAG],
                                 use_flair=config.props[ConfigKey.FLAIR_FLAG])
    dwi_processing = DWIProcessing()
    coregistration = Coregistration()
    tracts_generation = TractsGeneration()

    job_t1, job_aparc_aseg = t1_processing.add_t1_processing_steps(dax)
    job_b0, job_mask = dwi_processing.add_dwi_processing_steps(dax)
    job_t1_in_d, job_aparc_aseg_in_d = coregistration.add_coregistration_steps(dax, job_b0, job_t1,
                                                                               job_aparc_aseg)
    tracts_generation.add_tracts_generation_steps(dax, job_t1_in_d, job_mask, job_aparc_aseg_in_d)

    f = open(daxfile, "w")
    dax.writeXML(f)
    f.close()
