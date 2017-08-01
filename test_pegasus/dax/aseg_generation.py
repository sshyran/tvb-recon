from Pegasus.DAX3 import File, Job, Link
from mappings import T1Files, AsegGenJobNames, AsegFiles, Inputs, T1JobNames


class AsegGeneration(object):
    def __init__(self, lh_labels, rh_labels):
        self.lh_labels = lh_labels
        self.rh_labels = rh_labels

    def add_aseg_generation_steps(self, dax, job_recon):
        lh_aseg = File(AsegFiles.LH_ASEG.value)
        lh_aseg_annot = File(AsegFiles.LH_ASEG_ANNOT.value)
        fs_lut = File(Inputs.FS_LUT.value)
        job5 = Job(AsegGenJobNames.ASEG_CONCAT.value)
        job5.addArguments(lh_aseg, lh_aseg_annot, self.lh_labels, fs_lut)
        job5.uses(fs_lut, link=Link.INPUT)
        job5.uses(lh_aseg, link=Link.OUTPUT, transfer=False, register=False)
        job5.uses(lh_aseg_annot, link=Link.OUTPUT, transfer=False, register=False)
        dax.addJob(job5)

        rh_aseg = File(AsegFiles.RH_ASEG.value)
        rh_aseg_annot = File(AsegFiles.RH_ASEG_ANNOT.value)
        fs_lut = File(Inputs.FS_LUT.value)
        job6 = Job(AsegGenJobNames.ASEG_CONCAT.value)
        job6.addArguments(rh_aseg, rh_aseg_annot, self.rh_labels, fs_lut)
        job6.uses(fs_lut, link=Link.INPUT)
        job6.uses(rh_aseg, link=Link.OUTPUT, transfer=False, register=False)
        job6.uses(rh_aseg_annot, link=Link.OUTPUT, transfer=False, register=False)
        dax.addJob(job6)

        lbl_list = map(int, self.lh_labels.strip('"').split() + self.rh_labels.strip('"').split())
        for aseg_label in lbl_list:
            aparc_aseg_mgz = File(T1Files.APARC_ASEG_MGZ.value)
            norm_mgz = File(T1Files.NORM_MGZ.value)
            aseg_mgz = File(AsegFiles.ASEG_MGZ.value % aseg_label)
            job1 = Job(AsegGenJobNames.MRI_PRETESS.value)
            job1.addArguments(aparc_aseg_mgz, str(aseg_label), norm_mgz, aseg_mgz)
            job1.uses(aparc_aseg_mgz, link=Link.INPUT)
            job1.uses(norm_mgz, link=Link.INPUT)
            job1.uses(aseg_mgz, link=Link.OUTPUT, transfer=False, register=False)
            dax.addJob(job1)

            dax.depends(job1, job_recon)

            aseg_not_smooth = File(AsegFiles.ASEG_NOT_SMOOTH.value % aseg_label)
            job2 = Job(AsegGenJobNames.MRI_TESSELLATE.value)
            job2.addArguments(aseg_mgz, str(aseg_label), aseg_not_smooth)
            job2.uses(aseg_mgz, link=Link.INPUT)
            job2.uses(aseg_not_smooth, link=Link.OUTPUT, transfer=False, register=False)
            dax.addJob(job2)

            dax.depends(job2, job1)

            aseg_not_smooth_main = File(AsegFiles.ASEG_NOT_SMOOTH_MAIN.value % aseg_label)
            job3 = Job(AsegGenJobNames.MRIS_EXTRACT.value)
            job3.addArguments(aseg_not_smooth, aseg_not_smooth_main)
            job3.uses(aseg_not_smooth, link=Link.INPUT)
            job3.uses(aseg_not_smooth_main, link=Link.OUTPUT, transfer=False, register=False)
            dax.addJob(job3)

            dax.depends(job3, job2)

            aseg_lbl = File(AsegFiles.ASEG_LBL.value % aseg_label)
            job4 = Job(AsegGenJobNames.MRIS_SMOOTH.value)
            job4.addArguments("-nw", aseg_not_smooth_main, aseg_lbl)
            job4.uses(aseg_not_smooth_main, link=Link.INPUT)
            job4.uses(aseg_lbl, link=Link.OUTPUT, transfer=False, register=False)
            dax.addJob(job4)

            dax.depends(job4, job3)

            if aseg_label in map(int, self.lh_labels.strip('"').split()):
                job5.uses(aseg_lbl, link=Link.INPUT)
            else:
                job6.uses(aseg_lbl, link=Link.INPUT)
            dax.depends(job5, job4)
            dax.depends(job6, job4)

        lh_centered_aseg = File(AsegFiles.LH_CENTERED_ASEG.value)
        job7 = Job(T1JobNames.MRIS_CONVERT.value)
        job7.addArguments("--to-scanner", lh_aseg, lh_centered_aseg)
        job7.uses(lh_aseg, link=Link.INPUT)
        job7.uses(lh_centered_aseg, link=Link.OUTPUT, transfer=False, register=False)
        dax.addJob(job7)

        dax.depends(job7, job5)

        rh_centered_aseg = File(AsegFiles.RH_CENTERED_ASEG.value)
        job8 = Job(T1JobNames.MRIS_CONVERT.value)
        job8.addArguments("--to-scanner", rh_aseg, rh_centered_aseg)
        job8.uses(rh_aseg, link=Link.INPUT)
        job8.uses(rh_centered_aseg, link=Link.OUTPUT, transfer=False, register=False)
        dax.addJob(job8)

        dax.depends(job8, job6)

        return job7, job8
