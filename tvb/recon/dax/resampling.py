from Pegasus.DAX3 import File, Link, Job
from tvb.recon.dax import AtlasSuffix
from tvb.recon.dax.mappings import ResamplingFiles, T1Files, ResamplingJobNames, T1JobNames
from tvb.recon.dax.qc_snapshots import QCSnapshots


class Resampling(object):
    def __init__(self, subject, trg_subject, decim_factor, atlas_suffixes=[AtlasSuffix.DEFAULT]):
        self.subject = subject
        self.trg_subject = trg_subject
        self.decim_factor = decim_factor
        self.qc_snapshots = QCSnapshots.get_instance()
        self.atlas_suffixes = atlas_suffixes

    def add_surface_resampling_steps(self, dax, job_recon, job_t1):
        t1_mgz = File(T1Files.T1_MGZ.value)

        lh_pial = File(T1Files.LH_PIAL.value)
        rh_pial = File(T1Files.RH_PIAL.value)
        pials = [lh_pial, rh_pial]

        lh_pial_resamp = File(ResamplingFiles.LH_PIAL_RESAMP.value % self.trg_subject)
        rh_pial_resamp = File(ResamplingFiles.RH_PIAL_RESAMP.value % self.trg_subject)
        pials_resamp = [lh_pial_resamp, rh_pial_resamp]

        lh_aparc_annots = []
        rh_aparc_annots = []
        lh_aparc_annots_resamp = []
        rh_aparc_annots_resamp = []
        for atlas_suffix in self.atlas_suffixes:
            lh_aparc_annots.append(File(T1Files.LH_APARC_ANNOT.value % atlas_suffix))
            rh_aparc_annots.append(File(T1Files.RH_APARC_ANNOT.value % atlas_suffix))
            lh_aparc_annots_resamp.append(File(
                ResamplingFiles.LH_APARC_ANNOT_RESAMP.value % (self.trg_subject, atlas_suffix)))
            rh_aparc_annots_resamp.append(File(
                ResamplingFiles.RH_APARC_ANNOT_RESAMP.value % (self.trg_subject, atlas_suffix)))
        aparcs = [lh_aparc_annots, rh_aparc_annots]
        aparcs_resamp = [lh_aparc_annots_resamp, rh_aparc_annots_resamp]

        jobs1 = []
        jobs2 = [[], []]
        for ih, (hemi, h_pial, h_pial_resamp, h_aparc, h_aparc_resamp) in \
            enumerate(zip(["lh", "rh"], pials, pials_resamp,aparcs, aparcs_resamp)):
            jobs1.append(Job(ResamplingJobNames.MRI_SURF2SURF.value, node_label="Resample %s pial " % hemi))
            jobs1[-1].addArguments("pial", "--srcsubject", self.subject, "--trgsubject", self.trg_subject,
                                   "--hemi", hemi, "--sval-xyz", "pial", "--tval", "pial-%s" % self.trg_subject,
                                   "--tval-xyz", t1_mgz)
            jobs1[-1].uses(t1_mgz, link=Link.INPUT)
            jobs1[-1].uses(h_pial, link=Link.INPUT)
            jobs1[-1].uses(h_pial_resamp, link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(jobs1[-1])
            dax.depends(jobs1[-1], job_recon)

            for iatlas, (atlas_suffix, aparc, aparc_resamp) in \
                    enumerate(zip(self.atlas_suffixes, h_aparc, h_aparc_resamp)):
                jobs2[ih].append(Job(ResamplingJobNames.MRI_SURF2SURF.value,
                                     node_label="Resample %s_aparc%s" % (hemi, atlas_suffix)))
                # Correct for empty suffix for default atlas or argument $1 will be missing
                if len(atlas_suffix) == 0:
                    atlas_suffix_ = "default"
                else:
                    atlas_suffix_ = atlas_suffix
                jobs2[ih][-1].addArguments(atlas_suffix_,
                                           "--srcsubject", self.subject, "--trgsubject", self.trg_subject,
                                           "--hemi", hemi, "--sval-annot", aparc, "--tval", aparc_resamp)
                jobs2[ih][-1].uses(aparc, link=Link.INPUT)
                jobs2[ih][-1].uses(aparc_resamp, link=Link.OUTPUT, transfer=True, register=True)
                dax.addJob(jobs2[ih][-1])
                dax.depends(jobs2[ih][-1], job_recon)

        lh_centered_pial = File(ResamplingFiles.LH_CENTERED_PIAL_RESAMP.value % self.trg_subject)
        job5 = Job(T1JobNames.MRIS_CONVERT.value, node_label="mris_convert_lh")
        job5.addArguments("--to-scanner", lh_pial_resamp, lh_centered_pial)
        job5.uses(lh_pial_resamp, link=Link.INPUT)
        job5.uses(lh_centered_pial, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job5)
        dax.depends(job5, jobs1[0])

        rh_centered_pial = File(ResamplingFiles.RH_CENTERED_PIAL_RESAMP.value % self.trg_subject)
        job6 = Job(T1JobNames.MRIS_CONVERT.value, node_label="mris_convert_rh")
        job6.addArguments("--to-scanner", rh_pial_resamp, rh_centered_pial)
        job6.uses(rh_pial_resamp, link=Link.INPUT)
        job6.uses(rh_centered_pial, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job6)
        dax.depends(job6, jobs1[1])

        t1_nii_gz = File(T1Files.T1_NII_GZ.value)
        self.qc_snapshots.add_vol_surf_snapshot_step(dax, [job_t1, job5, job6],
                                                     t1_nii_gz, [lh_centered_pial, rh_centered_pial],
                                                     "resampl_pialsurf_t1")

        for iatlas, (atlas_suffix, lh_aparc_annot_resamp, rh_aparc_annot_resamp) in \
                enumerate(zip(self.atlas_suffixes, lh_aparc_annots_resamp, rh_aparc_annots_resamp)):
            self.qc_snapshots.add_surf_annot_snapshot_step(dax, [job5, jobs2[0][iatlas]],
                                                           lh_centered_pial, lh_aparc_annot_resamp,
                                                           "resampl_lh_pial%s" % atlas_suffix)
            self.qc_snapshots.add_surf_annot_snapshot_step(dax, [job6, jobs2[1][iatlas]],
                                                           rh_centered_pial, rh_aparc_annot_resamp,
                                                           "resampl_rh_pial%s" % atlas_suffix)

        jobs_pial = [job5, job6]
        jobs_aparc_per_atlas = [[] for _ in range(len(self.atlas_suffixes))]
        for jobs2h in jobs2:
            for iatlas, job2atlas in enumerate(jobs2h):
                jobs_aparc_per_atlas[iatlas].append(job2atlas)

        return jobs_pial, jobs_aparc_per_atlas

        # def add_annotation_resampling_steps(self, dax, job_aseg):
        # lh_aseg = File(AsegFiles.LH_ASEG.value)
        # rh_aseg = File(AsegFiles.RH_ASEG.value)
        # asegs = [lh_aseg, rh_aseg]

        # lh_aseg_annot = File(AsegFiles.LH_ASEG_ANNOT.value)
        # rh_aseg_annot = File(AsegFiles.RH_ASEG_ANNOT.value)
        # asegs = [lh_aseg_annot, rh_aseg_annot]

        # lh_aseg_resamp = File(ResamplingFiles.LH_ASEG_RESAMP.value % self.trg_subject)
        # rh_aseg_resamp = File(ResamplingFiles.RH_ASEG_RESAMP.value % self.trg_subject)
        # asegs_resamp = [lh_aseg_resamp, rh_aseg_resamp]

        # lh_aseg_annot_resamp = File(ResamplingFiles.LH_ASEG_ANNOT_RESAMP.value % self.trg_subject)
        # rh_aseg_annot_resamp = File(ResamplingFiles.RH_ASEG_ANNOT_RESAMP.value % self.trg_subject)
        # asegs_resamp = [lh_aseg_annot_resamp, rh_aseg_annot_resamp]

        # for idx, hemi in enumerate(["lh", "rh"]):
        # job2 = Job(ResamplingJobNames.MRIS_DECIMATE.value)
        # job2.addArguments("-d", self.decim_factor, asegs[idx], asegs_resamp[idx])
        # job2.uses(asegs[idx], link=Link.INPUT)
        # job2.uses(asegs_resamp[idx], link=Link.OUTPUT, transfer=True, register=True)
        # dax.addJob(job2)
        #
        # dax.depends(job2, job_aseg)

        # job2 = Job(ResamplingJobNames.?.value)
        # job2.addArguments()
        # job2.uses(asegs[idx], link=Link.INPUT)
        # job2.uses(asegs_resamp[idx], link=Link.OUTPUT, transfer=True, register=True)
        # dax.addJob(job2)
        #
        # dax.depends(job2, job_aseg)
