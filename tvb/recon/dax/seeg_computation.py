from Pegasus.DAX3 import File, Job, Link

from tvb.recon.dax.mappings import Inputs, SEEGCompFiles, T1JobNames, CoregJobNames, T1Files, SEEGCompJobNames


class SEEGComputation(object):
    def __init__(self, subj, ct_frmt, ct_elec_intensity_th):
        self.subj = subj
        self.ct_frmt = ct_frmt
        self.ct_elec_intensity_th = ct_elec_intensity_th

    def add_seeg_positions_computation_steps(self, dax):
        ct_input = File(Inputs.CT_INPUT.value)
        ct_ras = File(SEEGCompFiles.CT_RAS_NII_GZ.value)
        job1 = Job(T1JobNames.MRI_CONVERT.value, node_label="mri_convert CT to RAS")
        job1.addArguments(ct_input, ct_ras, "--out_orientation", "RAS")
        job1.uses(ct_input, Link.INPUT)
        job1.uses(ct_ras, Link.OUTPUT, register=True, transfer=True)
        dax.addJob(job1)

        t1_ras = File(T1Files.T1_NII_GZ.value)
        ct_in_t1 = File(SEEGCompFiles.CT_IN_T1_NII_GZ.value)
        ct_to_t1_mat = File(SEEGCompFiles.CT_TO_T1_MAT.value)
        job2 = Job(CoregJobNames.FLIRT.value, node_label="flirt register CT to T1")
        job2.addArguments(ct_ras, t1_ras, ct_to_t1_mat, ct_in_t1)
        job2.uses(t1_ras, Link.INPUT)
        job2.uses(ct_ras, Link.INPUT)
        job2.uses(ct_in_t1, Link.OUTPUT, transfer=True, register=True)
        job2.uses(ct_to_t1_mat, Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job2)

        dax.depends(job2, job1)

        brain_mgz = File(T1Files.BRAIN_MGZ.value)
        brain_ras = File(SEEGCompFiles.BRAIN_RAS_NII_GZ.value)
        job3 = Job(T1JobNames.MRI_CONVERT.value, node_label="mri_convert brain.mgz to RAS")
        job3.addArguments(brain_mgz, brain_ras, "--out_orientation", "RAS")
        job3.uses(brain_mgz, Link.INPUT)
        job3.uses(brain_ras, Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job3)

        brain_mask = File(SEEGCompFiles.BRAIN_MASK_NII_GZ.value)
        job4 = Job(SEEGCompJobNames.MRI_BINARIZE.value, node_label="mri_binarize brain.mgz to generate mask")
        job4.addArguments("--i", brain_ras, "--o", brain_mask, "--min", "10", "--erode", "4")
        job4.uses(brain_ras, Link.INPUT)
        job4.uses(brain_mask, Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job4)

        dax.depends(job4, job3)

        masked_ct = File(SEEGCompFiles.MASKED_CT_NII_GZ.value)
        job5 = Job(SEEGCompJobNames.MRI_BINARIZE.value, node_label="mri_binarize to mask CT with brain mask")
        job5.addArguments("--i", ct_in_t1, "--o", masked_ct, "--min", "1000", "--mask", brain_mask)
        job5.uses(ct_in_t1, Link.INPUT)
        job5.uses(brain_mask, Link.INPUT)
        job5.uses(masked_ct, Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job5)

        dax.depends(job5, job2)
        dax.depends(job5, job4)

        dilated_ct = File(SEEGCompFiles.DILATED_CT_NII_GZ.value)
        job6 = Job(SEEGCompJobNames.MRI_BINARIZE.value, node_label="mri_binarize to dilate & erode masked CT")
        job6.addArguments("--i", masked_ct, "--o", dilated_ct, "--min", "0.5", "--dilate", "2", "--erode", "1")
        job6.uses(masked_ct, Link.INPUT)
        job6.uses(dilated_ct, Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job6)

        dax.depends(job6, job5)

        labeled_ct = File(SEEGCompFiles.LABELED_CT_NII_GZ.value)
        job7 = Job(SEEGCompJobNames.LABEL_CT_WITH_DILATION.value, node_label="Label CT")
        job7.addArguments(masked_ct, dilated_ct, labeled_ct, self.subj)
        job7.uses(masked_ct, Link.INPUT)
        job7.uses(dilated_ct, Link.INPUT)
        job7.uses(labeled_ct, Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job7)

        dax.depends(job7, job6)

        schema_txt = File(Inputs.SCHEMA_TXT.value)
        job8 = Job(SEEGCompJobNames.GEN_SCHEMA_TXT.value, node_label="Generate seeg coords schema.txt")
        job8.addArguments(labeled_ct, schema_txt)
        job8.uses(labeled_ct, link=Link.INPUT)
        job8.uses(schema_txt, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job8)

        dax.depends(job8, job7)

        seeg_xyz = File(SEEGCompFiles.SEEG_XYZ.value)
        job9 = Job(SEEGCompJobNames.GEN_SEEG_XYZ.value, node_label="Generate seeg_xyz.txt")
        job9.addArguments(labeled_ct, schema_txt, seeg_xyz, self.subj)
        job9.uses(labeled_ct, Link.INPUT)
        job9.uses(schema_txt, Link.INPUT)
        job9.uses(seeg_xyz, Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job9)

        dax.depends(job9, job7)
        dax.depends(job9, job8)

        return job9