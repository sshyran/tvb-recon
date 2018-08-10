from Pegasus.DAX3 import File, Job, Link
from tvb.recon.dax.mappings import Inputs, DWIFiles, DWIJobNames
from tvb.recon.dax.qc_snapshots import QCSnapshots


class DWIProcessing(object):
    def __init__(self, dwi_reversed=False, dwi_frmt="mif", use_gradient=False, mrtrix_thrds="2", dwi_pe_dir="ap",
                 os="LINUX"):
        self.dwi_reversed = dwi_reversed
        self.dwi_format = dwi_frmt
        self.add_gradient = use_gradient
        self.mrtrix_threads = mrtrix_thrds
        self.dwi_pe_dir = dwi_pe_dir
        self.os = os
        self.qc_snapshots = QCSnapshots.get_instance()

    def add_dwi_processing_steps(self, dax):

        if self.dwi_reversed == "True":
            job1 = None
            job2 = None

            if self.dwi_format != "mif":

                if self.dwi_format == "dicom":
                    # TODO: is mrconvert interactive for reversed aquisition data? Should we use the next lines?
                    # mrchoose 0 mrconvert $DATA/DWI ./dwi_raw.mif -force
                    # mrchoose 1 mrconvert $DATA/DWI ./dwi_re_raw.mif -force
                    print("Not implemented!")

                else:
                    # TODO: to handle a possible case where gradient information is missing!
                    # Probably gradient (i.e., bvec, bval) is different for each direction...
                    dwi_input = File(DWIFiles.DWI_INPUT.value)
                    dwi_raw_mif = File(DWIFiles.DWI_RAW_MIF.value)
                    job1 = Job(DWIJobNames.MRCONVERT.value, node_label="mrconvert dwi input to dwi_raw.mif")
                    job1.addArguments(dwi_input)
                    job1.uses(dwi_input, link=Link.INPUT)
                    job1.uses(dwi_raw_mif, link=Link.OUTPUT, transfer=False, register=False)
                    dax.addJob(job1)

                    dwi_re_input = File(DWIFiles.DWI_RE_INPUT.value)
                    dwi_re_raw_mif = File(DWIFiles.DWI_RE_RAW_MIF.value)
                    job1re = Job(DWIJobNames.MRCONVERT.value,
                                  node_label="mrconvert dwi reversed input to dwi_re_raw.mif")
                    job1re.addArguments(dwi_re_input, dwi_re_raw_mif)
                    job1re.uses(dwi_re_input, link=Link.INPUT)
                    job1re.uses(dwi_re_raw_mif, link=Link.OUTPUT, transfer=True, register=False)
                    dax.addJob(job1re)

            dwi_pre_output = File(DWIFiles.DWI_MIF.value)
            job2 = Job(DWIJobNames.DWIPREPROC.value, node_label="DWI preprocessing (incl. reversed)")

            if self.os == "LINUX":
                job2.addArguments(dwi_raw_mif, dwi_pre_output, "-pe_dir", self.dwi_pe_dir, "-rpe_pair",
                                  dwi_raw_mif, dwi_re_raw_mif, "-nthreads", self.mrtrix_threads)
            else:
                job2.addArguments(self.dwi_pe_dir, dwi_raw_mif, dwi_pre_output, "-rpe_pair", dwi_raw_mif,
                                  dwi_re_raw_mif, "-nthreads",
                                  self.mrtrix_threads)
            job2.uses(dwi_raw_mif, link=Link.INPUT)
            job2.uses(dwi_raw_mif, link=Link.INPUT)
            job2.uses(dwi_pre_output, link=Link.OUTPUT, transfer=False, register=False)
            dax.addJob(job2)

            dax.depends(job2, job1)
            dax.depends(job2, job1re)

        else:

            dwi_input = File(Inputs.DWI_INPUT.value)
            dwi_raw_mif = File(DWIFiles.DWI_RAW_MIF.value)

            if self.add_gradient == "True":
                bvec_input = File(Inputs.DWI_BVEC.value)
                bval_input = File(Inputs.DWI_BVAL.value)

                job1 = Job(DWIJobNames.MRCONVERT.value,
                           node_label="mrconvert dwi input to dwi_raw.mif adding gradient bval & bvec")
                job1.addArguments(dwi_input, "-fsl", bvec_input, bval_input, dwi_input)
                job1.uses(dwi_input, link=Link.INPUT)
                job1.uses(bvec_input, link=Link.INPUT)
                job1.uses(bval_input, link=Link.INPUT)
                job1.uses(dwi_raw_mif, link=Link.OUTPUT, transfer=True, register=True)

            else:
                job1 = Job(DWIJobNames.MRCONVERT.value, node_label="mrconvert dwi input to dwi_raw.mif")
                job1.addArguments(dwi_input)
                job1.uses(dwi_input, link=Link.INPUT)
                job1.uses(dwi_raw_mif, link=Link.OUTPUT, transfer=True, register=True)

            dax.addJob(job1)

            dwi_pre_output = File(DWIFiles.DWI_MIF.value)
            job2 = Job(DWIJobNames.DWIPREPROC.value, node_label="DWI preprocessing")

            if self.os == "LINUX":
                job2.addArguments(dwi_raw_mif, dwi_pre_output, "-pe_dir", self.dwi_pe_dir, "-rpe_none", "-nthreads",
                                  self.mrtrix_threads)

            else:
                job2.addArguments(self.dwi_pe_dir, dwi_raw_mif, dwi_pre_output, "-rpe_none", "-nthreads",
                                  self.mrtrix_threads)

            job2.uses(dwi_raw_mif, link=Link.INPUT)
            job2.uses(dwi_pre_output, link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(job2)

            dax.depends(job2, job1)

        mask_output = File(DWIFiles.MASK_MIF.value)
        job3 = Job(DWIJobNames.DWI2MASK.value, node_label="Compute DWI mask")
        job3.addArguments(dwi_pre_output, mask_output, "-nthreads", self.mrtrix_threads)
        job3.uses(dwi_pre_output, link=Link.INPUT)
        job3.uses(mask_output, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job3)

        dax.depends(job3, job2)

        b0_output = File(DWIFiles.B0_NII_GZ.value)
        job4 = Job(DWIJobNames.DWIEXTRACT.value, node_label="Extract DWI B0")
        job4.addArguments(dwi_pre_output, b0_output)
        job4.uses(dwi_pre_output, link=Link.INPUT)
        job4.uses(b0_output, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job4)

        dax.depends(job4, job2)

        file_mask_nii_gz = File(DWIFiles.MASK_NII_GZ.value)
        job_convert_mask = Job(DWIJobNames.MRCONVERT.value, node_label="mrconvert mask mif->nii.gz")
        job_convert_mask.addArguments(mask_output, file_mask_nii_gz)
        job_convert_mask.uses(mask_output, link=Link.INPUT)
        job_convert_mask.uses(file_mask_nii_gz, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job_convert_mask)

        dax.depends(job_convert_mask, job3)

        self.qc_snapshots.add_2vols_snapshot_step(dax, [job_convert_mask, job4],
                                                  file_mask_nii_gz, b0_output, "mask_in_b0")

        return job4, job3
