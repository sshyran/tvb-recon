from Pegasus.DAX3 import File, Job, Link
from tvb.recon.dax import AtlasSuffix
from tvb.recon.dax.mappings import AsegFiles, AsegGenJobNames, Inputs, T1Files, T1JobNames, ResamplingFiles


class MappingDetails(object):
    def __init__(self, trg_subj, atlas_suffix=AtlasSuffix.DEFAULT):
        self.trg_subj = trg_subj
        self.atlas_suffix = atlas_suffix

    def add_mapping_details_step(self, dax, job_t1, job_lh_aseg, job_rh_aseg, jobs_resampling=None):

        # Input files:
        fs_lut = File(Inputs.FS_LUT.value)
        t1_file = File(T1Files.T1_NII_GZ.value)
        lh_aseg_annot = File(AsegFiles.LH_ASEG_ANNOT.value)
        rh_aseg_annot = File(AsegFiles.RH_ASEG_ANNOT.value)
        lh_subcort = File(AsegFiles.LH_CENTERED_ASEG.value)
        rh_subcort = File(AsegFiles.RH_CENTERED_ASEG.value)

        if jobs_resampling is None:
            lh_cortical_file = File(T1Files.LH_CENTERED_PIAL.value)
            rh_cortical_file = File(T1Files.RH_CENTERED_PIAL.value)
            lh_cortical_annot_file = File(T1Files.LH_APARC_ANNOT.value % self.atlas_suffix)
            rh_cortical_annot_file = File(T1Files.RH_APARC_ANNOT.value % self.atlas_suffix)
        else:
            lh_cortical_file = File(ResamplingFiles.LH_CENTERED_PIAL_RESAMP.value % self.trg_subj)
            rh_cortical_file = File(ResamplingFiles.RH_CENTERED_PIAL_RESAMP.value % self.trg_subj)
            lh_cortical_annot_file = File(
                ResamplingFiles.LH_APARC_ANNOT_RESAMP.value % (self.trg_subj, self.atlas_suffix))
            rh_cortical_annot_file = File(
                ResamplingFiles.RH_APARC_ANNOT_RESAMP.value % (self.trg_subj, self.atlas_suffix))

        # Output files:
        fs_custom = File(AsegFiles.FS_CUSTOM_TXT.value % self.atlas_suffix)
        centers = File(AsegFiles.CENTERS_TXT.value % self.atlas_suffix)
        areas = File(AsegFiles.AREAS_TXT.value % self.atlas_suffix)
        orientations = File(AsegFiles.ORIENTATIONS_TXT.value % self.atlas_suffix)
        cortical = File(AsegFiles.CORTICAL_TXT.value % self.atlas_suffix)
        rm_to_aparc_aseg = File(AsegFiles.RM_TO_APARC_ASEG_TXT.value % self.atlas_suffix)
        cort_region_mapping = File(AsegFiles.RM_CORT_TXT.value % self.atlas_suffix)
        subcort_region_mapping = File(AsegFiles.RM_SUBCORT_TXT.value % self.atlas_suffix)
        cort_surface = File(AsegFiles.SURF_CORT_ZIP.value)
        subcort_surface = File(AsegFiles.SURF_SUBCORT_ZIP.value)
        lh_dipoles = File(AsegFiles.LH_DIPOLES_TXT.value % self.atlas_suffix)
        rh_dipoles = File(AsegFiles.RH_DIPOLES_TXT.value % self.atlas_suffix)

        # Job config:
        job9 = Job(AsegGenJobNames.GEN_MAPPING_DETAILS.value)
        job9.addArguments(self.atlas_suffix, fs_lut, t1_file, lh_cortical_file, rh_cortical_file,
                          lh_cortical_annot_file, rh_cortical_annot_file, lh_subcort, rh_subcort, lh_aseg_annot,
                          rh_aseg_annot)
        job9.uses(fs_lut, link=Link.INPUT)
        job9.uses(t1_file, link=Link.INPUT)
        job9.uses(lh_cortical_annot_file, link=Link.INPUT)
        job9.uses(rh_cortical_annot_file, link=Link.INPUT)
        job9.uses(lh_cortical_file, link=Link.INPUT)
        job9.uses(rh_cortical_file, link=Link.INPUT)
        job9.uses(lh_aseg_annot, link=Link.INPUT)
        job9.uses(rh_aseg_annot, link=Link.INPUT)
        job9.uses(lh_subcort, link=Link.INPUT)
        job9.uses(rh_subcort, link=Link.INPUT)

        job9.uses(fs_custom, link=Link.OUTPUT, transfer=True, register=True)
        job9.uses(centers, link=Link.OUTPUT, transfer=True, register=True)
        job9.uses(areas, link=Link.OUTPUT, transfer=True, register=True)
        job9.uses(orientations, link=Link.OUTPUT, transfer=True, register=True)
        job9.uses(cortical, link=Link.OUTPUT, transfer=True, register=True)
        job9.uses(rm_to_aparc_aseg, link=Link.OUTPUT, transfer=True, register=True)
        job9.uses(cort_region_mapping, link=Link.OUTPUT, transfer=True, register=True)
        job9.uses(subcort_region_mapping, link=Link.OUTPUT, transfer=True, register=True)
        job9.uses(cort_surface, link=Link.OUTPUT, transfer=True, register=True)
        job9.uses(subcort_surface, link=Link.OUTPUT, transfer=True, register=True)
        job9.uses(lh_dipoles, link=Link.OUTPUT, transfer=True, register=True)
        job9.uses(rh_dipoles, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job9)
        dax.depends(job9, job_t1)
        dax.depends(job9, job_lh_aseg)
        dax.depends(job9, job_rh_aseg)

        if jobs_resampling is not None:
            for job_resampling in jobs_resampling:
                dax.depends(job9, job_resampling)

        return job9
