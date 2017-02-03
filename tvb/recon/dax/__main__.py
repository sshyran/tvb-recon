#!/usr/bin/env python

import time
import argparse
from Pegasus.DAX3 import ADAG
from tvb.recon.dax.config import Configuration
from tvb.recon.dax.transform import step_coregister_t1_dwi
from tvb.recon.dax.t1 import steps_recon_all
from tvb.recon.dax.diffusion import steps_dwi_preproc
from tvb.recon.dax.utils import write_dax


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a BNM flow")
    parser.add_argument("patient_file")
    args = parser.parse_args()

    dax = ADAG("BNM")
    dax.metadata("name", "Brain Network Model Reconstruction WorkFlow")
    dax.metadata("created-at", time.ctime())
    dax.metadata("flow-configuration", args.patient_file)
    config = Configuration(args.patient_file)

    relevant_t1_job = steps_recon_all(dax, config)
    relevant_dwi_job = steps_dwi_preproc(dax, config.diffusion)
    step_coregister_t1_dwi(dax, config, relevant_t1_job, relevant_dwi_job)

    write_dax(dax, config.main_dax_path)
