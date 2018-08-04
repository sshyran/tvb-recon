from Pegasus.DAX3 import File, Job, Link
from tvb.recon.dax import AtlasSuffix
from tvb.recon.dax.mappings import TractsGenFiles, OutputConvFiles, T1Files, AsegFiles


class OutputConversion(object):
    def __init__(self, atlas_suffixes=[AtlasSuffix.DEFAULT]):
        self.atlas_suffixes = atlas_suffixes

    def add_conversion_steps(self, dax, jobs_aparc_aseg, jobs_mapping_details, jobs_weights, jobs_lengths):

        weights_csv = []
        lenghts_csv = []
        centers = []
        areas = []
        orientations = []
        cortical = []
        rm_to_aparc_aseg = []
        jobs = []
        for atlas_suffix, job_aparc_aseg, job_mapping_details, job_weights, job_lengths in \
            zip(self.atlas_suffixes, jobs_aparc_aseg, jobs_mapping_details, jobs_weights, jobs_lengths):

            weights_csv.append(File(TractsGenFiles.TRACT_COUNTS.value % atlas_suffix))
            lenghts_csv.append(File(TractsGenFiles.TRACT_LENGHTS.value % atlas_suffix))

            centers.append(File(AsegFiles.CENTERS_TXT.value % atlas_suffix))
            areas.append(File(AsegFiles.AREAS_TXT.value % atlas_suffix))
            orientations.append(File(AsegFiles.ORIENTATIONS_TXT.value % atlas_suffix))
            cortical.append(File(AsegFiles.CORTICAL_TXT.value % atlas_suffix))
            rm_to_aparc_aseg.append(File(AsegFiles.RM_TO_APARC_ASEG_TXT.value %atlas_suffix))

            # aparc_aseg = File(T1Files.APARC_ASEG_NII_GZ.value)

            if len(atlas_suffix) == 0:
                atlas_name = "default"
            else:
                atlas_name = atlas_suffix[1:]
            jobs.append(Job("convert_output", node_label="convert_output of atlas " + atlas_name))
            jobs[-1].addArguments(weights_csv[-1], lenghts_csv[-1], atlas_suffix)

            jobs[-1].uses(weights_csv[-1], link=Link.INPUT)
            jobs[-1].uses(lenghts_csv[-1], link=Link.INPUT)
            jobs[-1].uses(centers[-1], link=Link.INPUT)
            jobs[-1].uses(areas[-1], link=Link.INPUT)
            jobs[-1].uses(orientations[-1], link=Link.INPUT)
            jobs[-1].uses(cortical[-1], link=Link.INPUT)
            jobs[-1].uses(rm_to_aparc_aseg[-1], link=Link.INPUT)
            # jobs[-1].uses(aparc_aseg, link=Link.INPUT)

            jobs[-1].uses(File(OutputConvFiles.APARC_ASEG_COR_NII_GZ.value % atlas_suffix),
                          link=Link.OUTPUT, transfer=True, register=False)
            jobs[-1].uses(File(OutputConvFiles.CONNECTIVITY_ZIP.value % atlas_suffix),
                          link=Link.OUTPUT, transfer=True, register=False)

            jobs[-1].uses(File(T1Files.T1_NII_GZ.value), link=Link.INPUT)
            jobs[-1].uses(File(T1Files.APARC_ASEG_NII_GZ.value % atlas_suffix), link=Link.INPUT)

            dax.addJob(jobs[-1])

            dax.depends(jobs[-1], job_aparc_aseg)
            dax.depends(jobs[-1], job_mapping_details)
            dax.depends(jobs[-1], job_weights)
            dax.depends(jobs[-1], job_lengths)

        return jobs
