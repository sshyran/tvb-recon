from Pegasus.DAX3 import File, Job, Link
from mappings import TractsGenFiles, Inputs, T1Files, AsegFiles, OutputConvFiles


class OutputConversion(object):
    def add_conversion_steps(self, dax, job_aparc_aseg, job_aseg_lh, job_aseg_rh):
        weights_csv = File(TractsGenFiles.TRACT_COUNTS.value)
        lenghts_csv = File(TractsGenFiles.TRACT_LENGHTS.value)
        fs_lut = File(Inputs.FS_LUT.value)

        job26 = Job("convert_output")
        job26.addArguments(weights_csv, lenghts_csv, fs_lut)

        job26.uses(weights_csv, link=Link.INPUT)
        job26.uses(lenghts_csv, link=Link.INPUT)
        job26.uses(fs_lut, link=Link.INPUT)

        job26.uses(File(OutputConvFiles.RM_CORT_TXT.value), link=Link.OUTPUT, transfer=True, register=False)
        job26.uses(File(OutputConvFiles.RM_SUBCORT_TXT.value), link=Link.OUTPUT, transfer=True, register=False)
        job26.uses(File(OutputConvFiles.SURF_CORT_ZIP.value), link=Link.OUTPUT, transfer=True, register=False)
        job26.uses(File(OutputConvFiles.SURF_SUBCORT_ZIP.value), link=Link.OUTPUT, transfer=True, register=False)
        job26.uses(File(OutputConvFiles.APARC_ASEG_COR_NII_GZ.value), link=Link.OUTPUT, transfer=True, register=False)
        job26.uses(File(OutputConvFiles.CONNECTIVITY_ZIP.value), link=Link.OUTPUT, transfer=True, register=False)

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
