from Pegasus.DAX3 import File, Job, Link
from tvb.recon.dax.mappings import CoregFiles, TractsGenFiles, TractsGenJobNames, DWIJobNames, DWIFiles
from tvb.recon.dax.qc_snapshots import QCSnapshots


class TractsGeneration(object):
    def __init__(self, dwi_multi_shell=False, mrtrix_threads="2", strmlns_no="25M", strmlns_sift_no="5M",
                 strmlns_size="250", strmlns_step="0.5", os="LINUX"):
        self.dwi_multi_shell = dwi_multi_shell
        self.mrtrix_threads = mrtrix_threads
        self.strmlns_no = strmlns_no
        self.strmlns_sift_no = strmlns_sift_no
        self.strmlns_size = strmlns_size
        self.strmlns_step = strmlns_step
        self.qc_snapshots = QCSnapshots.get_instance()
        self.os = os

    # job_t1_in_d = job12, job_mask = job5, job_aparc_aseg_in_d = job21
    def add_tracts_generation_steps(self, dax, job_t1_in_d, job_mask):

        #-------------------------------------------Tractography--------------------------------------------------------

        t1_in_d = File(CoregFiles.T1_IN_D.value)
        file_5tt = File(TractsGenFiles.FILE_5TT_MIF.value)
        job1 = Job(TractsGenJobNames.JOB_5TTGEN.value, node_label="Generate 5tt MIF")
        job1.addArguments(t1_in_d, file_5tt)
        job1.uses(t1_in_d, link=Link.INPUT)
        job1.uses(file_5tt, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job1)

        dax.depends(job1, job_t1_in_d)

        file_gmwmi = File(TractsGenFiles.GMWMI_MIF.value)
        job2 = Job(TractsGenJobNames.JOB_5TT2GMWMI.value, node_label="Extract GMWMI")
        job2.addArguments(file_5tt, file_gmwmi, "-nthreads", self.mrtrix_threads)
        job2.uses(file_5tt, link=Link.INPUT)
        job2.uses(file_gmwmi, link=Link.OUTPUT, transfer=False, register=False)
        dax.addJob(job2)

        dax.depends(job2, job1)

        file_gmwmi_nii_gz = File(TractsGenFiles.GMWMI_NII_GZ.value)
        job_gmwmi_convert = Job(DWIJobNames.MRCONVERT.value, node_label="mrconvert gmwmi mif->nii")
        job_gmwmi_convert.addArguments(file_gmwmi, file_gmwmi_nii_gz)
        job_gmwmi_convert.uses(file_gmwmi, link=Link.INPUT)
        job_gmwmi_convert.uses(file_gmwmi_nii_gz, link=Link.OUTPUT, transfer=False, register=False)
        dax.addJob(job_gmwmi_convert)

        dax.depends(job_gmwmi_convert, job2)

        self.qc_snapshots.add_2vols_snapshot_step(dax, [job_gmwmi_convert],
                                                  t1_in_d, file_gmwmi_nii_gz, "gmwmi_in_t1_in_d")

        file_5ttvis = File(TractsGenFiles.FILE_5TTVIS_MIF.value)
        job3 = Job(TractsGenJobNames.JOB_5TT2VIS.value, node_label="Generate TT2VIS MIF")
        job3.addArguments(file_5tt, file_5ttvis)
        job3.uses(file_5tt, link=Link.INPUT)
        job3.uses(file_5ttvis, link=Link.OUTPUT, transfer=False, register=False)
        dax.addJob(job3)

        dax.depends(job3, job2)

        dwi_mif = File(DWIFiles.DWI_MIF.value)
        mask_mif = File(DWIFiles.MASK_MIF.value)

        if self.dwi_multi_shell == "True":
            file_RF_WM = File(TractsGenFiles.RF_WM.value)
            file_RF_GM = File(TractsGenFiles.RF_GM.value)
            file_RF_CSF = File(TractsGenFiles.RF_CSF.value)
            file_RF_voxels = File(TractsGenFiles.RF_VOXELS.value)
            job4 = Job(TractsGenJobNames.DWI2RESPONSE_MSMT.value, node_label="Compute MSMT DWI response")
            job4.addArguments(dwi_mif, file_5tt, file_RF_WM, file_RF_GM, file_RF_CSF,
                              file_RF_voxels, self.mrtrix_threads)
            job4.uses(dwi_mif, link=Link.INPUT)
            job4.uses(file_5tt, link=Link.INPUT)
            job4.uses(file_RF_WM, link=Link.OUTPUT, transfer=True, register=False)
            job4.uses(file_RF_GM, link=Link.OUTPUT, transfer=True, register=False)
            job4.uses(file_RF_CSF, link=Link.OUTPUT, transfer=True, register=False)
            job4.uses(file_RF_voxels, link=Link.OUTPUT, transfer=True, register=False)
            dax.addJob(job4)

            dax.depends(job4, job3)

            gm_mif = File(DWIFiles.GM_MIF.value)
            csf_mif = File(DWIFiles.CSF_MIF.value)
            file_wm_fod = File(TractsGenFiles.WM_FOD_MIF.value)
            # TODO: does msdwi2fod exist? should we use dwi2fod with the same args?
            job5 = Job(TractsGenJobNames.MSDWI2FOD.value, node_label="Compute MSMT WM FOD")
            job5.addArguments("msmt_csd", dwi_mif, file_RF_WM, file_wm_fod, file_RF_GM, gm_mif, file_RF_CSF, csf_mif,
                              "-mask", mask_mif, "-nthreads", self.mrtrix_threads)
            job5.uses(dwi_mif, link=Link.INPUT)
            job5.uses(file_RF_WM, link=Link.INPUT)
            job5.uses(file_RF_GM, link=Link.INPUT)
            job5.uses(file_RF_CSF, link=Link.INPUT)
            job5.uses(mask_mif, link=Link.INPUT)
            job5.uses(file_wm_fod, link=Link.OUTPUT, transfer=True, register=False)
            job5.uses(gm_mif, link=Link.OUTPUT, transfer=True, register=False)
            job5.uses(csf_mif, link=Link.OUTPUT, transfer=True, register=False)
            dax.addJob(job5)

            dax.depends(job5, job4)

            last_job = job5

        else:
            file_response = File(TractsGenFiles.RESPONSE_TXT.value)
            job4 = Job(TractsGenJobNames.DWI2RESPONSE.value, node_label="Compute DWI response")
            job4.addArguments(dwi_mif, file_response, mask_mif)
            job4.uses(dwi_mif, link=Link.INPUT)
            job4.uses(mask_mif, link=Link.INPUT)
            job4.uses(file_response, link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(job4)

            dax.depends(job4, job_mask)

            file_wm_fod = File(TractsGenFiles.WM_FOD_MIF.value)
            job5 = Job(TractsGenJobNames.DWI2FOD.value, node_label="Compute WM FOD")
            job5.addArguments("csd", dwi_mif, file_response, file_wm_fod, "-mask", mask_mif,
                              "-nthreads",
                              self.mrtrix_threads)
            job5.uses(dwi_mif, link=Link.INPUT)
            job5.uses(file_response, link=Link.INPUT)
            job5.uses(mask_mif, link=Link.INPUT)
            job5.uses(file_wm_fod, link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(job5)

            dax.depends(job5, job4)

            last_job = job5

        file_strmlns = File(TractsGenFiles.FILE_TCK.value % self.strmlns_no)
        job6 = Job(TractsGenJobNames.TCKGEN.value, node_label="Generate tracts")

        if self.os == "LINUX":
            job6.addArguments(file_wm_fod, file_strmlns, "-select", self.strmlns_no, "-seed_gmwmi", file_gmwmi, "-act",
                              file_5tt, "-seed_unidirectional", "-maxlength", self.strmlns_size, "-step",
                              self.strmlns_step, "-nthreads", self.mrtrix_threads)
        else:
            job6.addArguments(file_wm_fod, file_strmlns, "-number", self.strmlns_no, "-seed_gmwmi", file_gmwmi, "-act",
                              file_5tt, "-unidirectional", "-maxlength", self.strmlns_size, "-step", self.strmlns_step,
                              "-nthreads", self.mrtrix_threads)
        job6.uses(file_wm_fod, link=Link.INPUT)
        job6.uses(file_gmwmi, link=Link.INPUT)
        job6.uses(file_5tt, link=Link.INPUT)
        job6.uses(file_strmlns, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job6)

        dax.depends(job6, last_job)
        dax.depends(job6, job1)

        file_strmlns_sift = File(TractsGenFiles.FILE_SIFT_TCK.value % self.strmlns_sift_no)
        job_tcksift = Job(TractsGenJobNames.TCKSIFT.value, node_label="Tracts SIFT")
        job_tcksift.addArguments(file_strmlns, file_wm_fod, file_strmlns_sift, "-term_number", self.strmlns_sift_no, "-act",
                          file_5tt, "-nthreads", self.mrtrix_threads)
        job_tcksift.uses(file_strmlns, link=Link.INPUT)
        job_tcksift.uses(file_wm_fod, link=Link.INPUT)
        job_tcksift.uses(file_5tt, link=Link.INPUT)
        job_tcksift.uses(file_strmlns_sift, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job_tcksift)

        dax.depends(job_tcksift, job6)
        dax.depends(job_tcksift, job1)

        # TODO: find out if job8 is necessary for parcellation-wise (i.e., non voxel-wise) connectomes

        b0_nii_gz = File(DWIFiles.B0_NII_GZ.value)
        file_tdi_ends = File(TractsGenFiles.TDI_ENDS_MIF.value)
        job8 = Job(TractsGenJobNames.TCKMAP.value, node_label="TCKMAP tracts")
        job8.addArguments(file_strmlns_sift, file_tdi_ends, "-vox", "1", "-template", b0_nii_gz)
        job8.uses(file_strmlns_sift, link=Link.INPUT)
        job8.uses(b0_nii_gz, link=Link.INPUT)
        job8.uses(file_tdi_ends, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job8)

        dax.depends(job8, job_tcksift)

        file_tdi_ends_nii_gz = File(TractsGenFiles.TDI_ENDS_NII_GZ.value)
        job_convert_tdi_ends = Job(DWIJobNames.MRCONVERT.value, node_label="mrconvert tdi_ends mif->nii")
        job_convert_tdi_ends.addArguments(file_tdi_ends, file_tdi_ends_nii_gz)
        job_convert_tdi_ends.uses(file_tdi_ends, link=Link.INPUT)
        job_convert_tdi_ends.uses(file_tdi_ends_nii_gz, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job_convert_tdi_ends)

        dax.depends(job_convert_tdi_ends, job8)

        self.qc_snapshots.add_2vols_snapshot_step(dax, [job_convert_tdi_ends],
                                                  t1_in_d, file_tdi_ends_nii_gz, "tdi_ends_in_t1_in_d")

        return job_tcksift
