import sys

import time
from Pegasus.DAX3 import ADAG, Job, Link, File
from t1_processing import T1Processing
from dwi_processing import DWIProcessing
from coregistration import Coregistration
from configuration import Configuration, ConfigKey
from aseg_generation import AsegGeneration
from tracts_generation import TractsGeneration
from mappings import TractsGenFiles, Inputs, T1Files, AsegFiles

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: %s DAXFILE\n" % (sys.argv[0]))
        sys.exit(1)
    daxfile = sys.argv[1]
    patient_file = sys.argv[2]

    dax = ADAG("TVB-PIPELINE")
    dax.metadata("created", time.ctime())

    config = Configuration(patient_file)

    t1_processing = T1Processing(t1_frmt=config.props[ConfigKey.T1_FRMT], use_t2=config.props[ConfigKey.T2_FLAG],
                                 use_flair=config.props[ConfigKey.FLAIR_FLAG])
    dwi_processing = DWIProcessing()
    coregistration = Coregistration()
    tracts_generation = TractsGeneration()
    aseg_generation = AsegGeneration(config.props[ConfigKey.ASEG_LH_LABELS], config.props[ConfigKey.ASEG_RH_LABELS])

    job_t1, job_aparc_aseg = t1_processing.add_t1_processing_steps(dax)
    job_b0, job_mask = dwi_processing.add_dwi_processing_steps(dax)
    job_t1_in_d, job_aparc_aseg_in_d = coregistration.add_coregistration_steps(dax, job_b0, job_t1,
                                                                               job_aparc_aseg)
    tracts_generation.add_tracts_generation_steps(dax, job_t1_in_d, job_mask, job_aparc_aseg_in_d)
    job_aseg_lh, job_aseg_rh = aseg_generation.add_aseg_generation_steps(dax, job_aparc_aseg)

    weights_csv = File(TractsGenFiles.TRACT_COUNTS.value)
    lenghts_csv = File(TractsGenFiles.TRACT_LENGHTS.value)
    fs_lut = File(Inputs.FS_LUT.value)
    job26 = Job("convert_output")
    job26.addArguments(weights_csv, lenghts_csv, fs_lut)
    job26.uses(weights_csv, link=Link.INPUT)
    job26.uses(lenghts_csv, link=Link.INPUT)
    job26.uses(fs_lut, link=Link.INPUT)
    job26.uses(File("region_mapping_cort.txt"), link=Link.OUTPUT, transfer=True, register=False)
    job26.uses(File("region_mapping_subcort.txt"), link=Link.OUTPUT, transfer=True, register=False)
    job26.uses(File("surface_cort.zip"), link=Link.OUTPUT, transfer=True, register=False)
    job26.uses(File("surface_subcort.zip"), link=Link.OUTPUT, transfer=True, register=False)
    job26.uses(File("aparc+aseg-cor.nii.gz"), link=Link.OUTPUT, transfer=True, register=False)
    job26.uses(File("connectivity.zip"), link=Link.OUTPUT, transfer=True, register=False)

    job26.uses(File(T1Files.T1_NII_GZ.value), link=Link.INPUT)
    job26.uses(File(T1Files.APARC_ASEG_NII_GZ.value), link=Link.INPUT)
    job26.uses(File(T1Files.LH_CENTERED_PIAL.value), link=Link.INPUT)
    job26.uses(File(T1Files.RH_CENTERED_PIAL.value), link=Link.INPUT)
    job26.uses(File(T1Files.LH_APARC_ANNOT.value), link=Link.INPUT)
    job26.uses(File(T1Files.RH_APARC_ANNOT.value), link=Link.INPUT)
    job26.uses(File(AsegFiles.LH_CENTERED_ASEG.value), link=Link.INPUT)
    job26.uses(File(AsegFiles.RH_CENTERED_ASEG.value), link=Link.INPUT)
    job26.uses(File(AsegFiles.LH_ASEG_ANNOT.value), link=Link.INPUT)
    job26.uses(File(AsegFiles.RH_ASEG_ANNOT.value), link=Link.INPUT)

    dax.addJob(job26)

    dax.depends(job26, job_aparc_aseg)
    dax.depends(job26, job_aseg_lh)
    dax.depends(job26, job_aseg_rh)

    f = open(daxfile, "w")
    dax.writeXML(f)
    f.close()
