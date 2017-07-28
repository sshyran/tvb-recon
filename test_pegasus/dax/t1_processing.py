from Pegasus.DAX3 import File, Job, Link
from mappings import Inputs, T1Files, T1JobNames


class T1Processing(object):
    def __init__(self, t1_frmt="nii", use_t2=False, t2_frmt="nii", use_flair=False, flair_frmt="nii", openmp_thrds="4"):
        self.t1_format = t1_frmt
        self.t2_flag = use_t2
        self.t2_format = t2_frmt
        self.flair_flag = use_flair
        self.flair_format = flair_frmt
        self.openmp_threads = openmp_thrds

    def _ensure_input_format(self, file_format, input_name, output_name, dax):
        input_file = File(input_name)
        output_file = None
        job = None

        if file_format == "dicom":
            output_file = File(output_name)
            job = Job(T1JobNames.MRI_CONVERT.value)
            job.addArguments(input_file, output_file)
            job.uses(input_file, link=Link.INPUT)
            job.uses(output_file, link=Link.OUTPUT, transfer=False, register=False)
            dax.addJob(job)

        if output_file is None:
            output_file = input_file

        return output_file, job

    def add_t1_processing_steps(self, dax):
        # Here we should take from configuration the:
        # SUBJ_NAME, T1,T2,FLAIR format, T2,FLAIR use, OPENMP threads

        subject = "TVB2PEG23"

        t1_input = Inputs.T1_INPUT.value
        t1_converted = T1Files.T1_INPUT_CONVERTED.value
        t1_output, job1 = self._ensure_input_format(self.t1_format, t1_input, t1_converted, dax)

        t1_mgz_output = File(T1Files.T1_MGZ.value)
        aparc_aseg_mgz_vol = File(T1Files.APARC_ASEG_MGZ.value)
        job2 = Job(T1JobNames.RECON_ALL.value, node_label="Recon-all for T1")
        job2.addArguments(subject, t1_output, self.openmp_threads)
        job2.uses(t1_output, link=Link.INPUT)
        job2.uses(t1_mgz_output, link=Link.OUTPUT, transfer=True, register=False)
        job2.uses(aparc_aseg_mgz_vol, link=Link.OUTPUT, transfer=True, register=False)
        dax.addJob(job2)

        if job1 is not None:
            dax.depends(job2, job1)

        last_job = job2

        if self.t2_flag is True:
            t2_in = Inputs.T2_INPUT.value
            t2_converted = T1Files.T2_CONVERTED.value
            t2_input, job_convert = self._ensure_input_format(self.t2_format, t2_in, t2_converted, dax)

            job = Job("autorecon3-t2")
            job.addArguments(subject, t2_input, self.openmp_threads)
            job.uses(t2_input, link=Link.INPUT)
            job.uses(aparc_aseg_mgz_vol, link=Link.OUTPUT, transfer=True, register=False)
            dax.addJob(job)

            if job_convert is not None:
                dax.depends(job, job_convert)
            dax.depends(job, last_job)

            last_job = job

        if self.flair_flag is True:
            flair_in = Inputs.FLAIR_INPUT.value
            flair_converted = T1Files.FLAIR_CONVERTED.value
            flair_input, job_convert = self._ensure_input_format(self.flair_format, flair_in, flair_converted, dax)

            job = Job("autorecon3-flair")
            job.addArguments(subject, flair_input, self.openmp_threads)
            job.uses(flair_input, link=Link.INPUT)
            job.uses(aparc_aseg_mgz_vol, link=Link.OUTPUT, transfer=True, register=False)
            dax.addJob(job)

            if job_convert is not None:
                dax.depends(job, job_convert)
            dax.depends(job, last_job)

            last_job = job

        t1_nii_gz_vol = File(T1Files.T1_NII_GZ.value)
        job3 = Job(T1JobNames.MRI_CONVERT.value, node_label="Convert T1 to NIFTI with good orientation")
        job3.addArguments(t1_mgz_output, t1_nii_gz_vol, "--out_orientation", "RAS")
        job3.uses(t1_mgz_output, link=Link.INPUT)
        job3.uses(t1_nii_gz_vol, link=Link.OUTPUT, transfer=True, register=False)
        dax.addJob(job3)

        dax.depends(job3, last_job)

        aparc_aseg_nii_gz_vol = File(T1Files.APARC_ASEG_NII_GZ.value)
        job4 = Job(T1JobNames.MRI_CONVERT.value, node_label="Convert APARC+ASEG to NIFTI with good orientation")
        job4.addArguments(aparc_aseg_mgz_vol, aparc_aseg_nii_gz_vol, "--out_orientation", "RAS", "-rt", "nearest")
        job4.uses(aparc_aseg_mgz_vol, link=Link.INPUT)
        job4.uses(aparc_aseg_nii_gz_vol, link=Link.OUTPUT, transfer=True, register=False)
        dax.addJob(job4)

        dax.depends(job4, last_job)

        return job3, job4
