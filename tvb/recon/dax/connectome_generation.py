from Pegasus.DAX3 import File, Job, Link
from tvb.recon.dax import AtlasSuffix
from tvb.recon.dax.mappings import CoregFiles, TractsGenFiles, TractsGenJobNames, AsegFiles, Inputs
from tvb.recon.dax.qc_snapshots import QCSnapshots


class ConnectomeGeneration(object):
    def __init__(self, strmlns_sift_no="5M", atlas_suffix=AtlasSuffix.DEFAULT):
        self.strmlns_sift_no = strmlns_sift_no
        self.atlas_suffix = atlas_suffix
        self.qc_snapshots = QCSnapshots.get_instance()

    # job_t1_in_d = job12, job_mask = job5, job_aparc_aseg_in_d = job21
    def add_tracts_generation_steps(self, dax, job_tcksift, job_aparc_aseg_in_d, job_fs_custom):

        t1_in_d = File(CoregFiles.T1_IN_D.value)
        file_strmlns_sift = File(TractsGenFiles.FILE_SIFT_TCK.value % self.strmlns_sift_no)

        fs_custom = File(AsegFiles.FS_CUSTOM_TXT.value % self.atlas_suffix)
        aparc_aseg_in_d = File(CoregFiles.APARC_ASEG_IN_D.value % self.atlas_suffix)
        file_vol_lbl = File(TractsGenFiles.VOLUME_LBL_NII_GZ.value % self.atlas_suffix)
        fs_color_lut = File(Inputs.FS_LUT.value)
        job9 = Job(TractsGenJobNames.LABEL_CONVERT.value,
                   node_label="Compute labeled APARC%s+ASEG for tracts" % self.atlas_suffix)
        job9.addArguments(aparc_aseg_in_d, fs_color_lut, fs_custom, file_vol_lbl)
        job9.uses(aparc_aseg_in_d, link=Link.INPUT)
        job9.uses(fs_color_lut, link=Link.INPUT)
        job9.uses(fs_custom, link=Link.INPUT)
        job9.uses(file_vol_lbl, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job9)

        dax.depends(job9, job_fs_custom)
        dax.depends(job9, job_aparc_aseg_in_d)

        self.qc_snapshots.add_2vols_snapshot_step(dax, [job9], t1_in_d, file_vol_lbl,
                                                  "labelled_aparc%s+aseg_in_t1_in_d" % self.atlas_suffix)

        if len(self.atlas_suffix) == 0:
            atlas_name = "default"
        else:
            atlas_name = self.atlas_suffix[1:]

        file_aparc_aseg_counts_csv = File(TractsGenFiles.TRACT_COUNTS.value % self.atlas_suffix)
        job10 = Job(TractsGenJobNames.TCK2CONNECTOME.value,
                    node_label="Generate weigths for atlas %s "% atlas_name)
        job10.addArguments(file_strmlns_sift, file_vol_lbl, "-assignment_radial_search", "2",
                                file_aparc_aseg_counts_csv)
        job10.uses(file_strmlns_sift, link=Link.INPUT)
        job10.uses(file_vol_lbl, link=Link.INPUT)
        job10.uses(file_aparc_aseg_counts_csv, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job10)

        dax.depends(job10, job_tcksift)
        dax.depends(job10, job9)

        file_aparc_aseg_mean_tract_lengths_csv = File(TractsGenFiles.TRACT_LENGHTS.value % self.atlas_suffix)
        job11 = Job(TractsGenJobNames.TCK2CONNECTOME.value,
                    node_label="Generate tract lengths for atlas %s" % atlas_name)
        job11.addArguments(file_strmlns_sift, file_vol_lbl,
                                "-assignment_radial_search", "2", "-scale_length",
                                "-stat_edge", "mean", file_aparc_aseg_mean_tract_lengths_csv)
        job11.uses(file_strmlns_sift, link=Link.INPUT)
        job11.uses(file_vol_lbl, link=Link.INPUT)
        job11.uses(file_aparc_aseg_mean_tract_lengths_csv, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job11)

        dax.depends(job11, job_tcksift)
        dax.depends(job11, job9)

        return job10, job11
