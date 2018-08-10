from Pegasus.DAX3 import Job, File, Link
from tvb.recon.dax import AtlasSuffix
from tvb.recon.dax.mappings import T1Files, ResamplingJobNames, SourceModelFiles, HeadModelJobNames, AsegFiles, \
    HeadModelFiles


class SourceModel(object):
    def __init__(self, subject, trg_subject, atlas_suffixes=[AtlasSuffix.DEFAULT]):
        self.subject = subject
        self.trg_subject = trg_subject
        self.atlas_suffixes = atlas_suffixes

    def add_source_model_steps(self, dax, job_head_model, jobs_mapping_details):
        t1_mgz = File(T1Files.T1_MGZ.value)
        lh_white = File(T1Files.LH_WHITE.value)
        rh_white = File(T1Files.RH_WHITE.value)
        whites = [lh_white, rh_white]

        head_model_geom = File(HeadModelFiles.HEAD_MODEL_GEOM.value)
        head_model_cond = File(HeadModelFiles.HEAD_MODEL_COND.value)

        bem_tri_surfs = [File(HeadModelFiles.INNER_SKULL_SURFACE_LOW_TRI.value % self.subject),
                         File(HeadModelFiles.OUTER_SKULL_SURFACE_LOW_TRI.value % self.subject),
                         File(HeadModelFiles.OUTER_SKIN_SURFACE_LOW_TRI.value % self.subject),
                         File(HeadModelFiles.BRAIN_SURFACE_LOW_TRI.value % self.subject)]

        lh_white_resamp = File(SourceModelFiles.LH_WHITE_RESAMP.value % self.trg_subject)
        rh_white_resamp = File(SourceModelFiles.RH_WHITE_RESAMP.value % self.trg_subject)
        whites_resamp = [lh_white_resamp, rh_white_resamp]

        lh_white_tri = File(SourceModelFiles.LH_WHITE_RESAMP_TRI.value % self.trg_subject)
        rh_white_tri = File(SourceModelFiles.RH_WHITE_RESAMP_TRI.value % self.trg_subject)
        whites_resamp_tri = [lh_white_tri, rh_white_tri]

        lh_white_ssm = File(SourceModelFiles.LH_WHITE_RESAMP_SSM.value % self.trg_subject)
        rh_white_ssm = File(SourceModelFiles.RH_WHITE_RESAMP_SSM.value % self.trg_subject)
        whites_resamp_ssm = [lh_white_ssm, rh_white_ssm]

        lh_dipoles_file = []
        rh_dipoles_file = []
        lh_white_dsm = []
        rh_white_dsm = []
        for atlas_suffix in self.atlas_suffixes:
            lh_dipoles_file.append(File(AsegFiles.LH_DIPOLES_TXT.value % atlas_suffix))
            rh_dipoles_file.append(File(AsegFiles.RH_DIPOLES_TXT.value % atlas_suffix))
            lh_white_dsm.append(File(SourceModelFiles.LH_WHITE_RESAMP_DSM.value % (self.trg_subject, atlas_suffix)))
            rh_white_dsm.append(File(SourceModelFiles.RH_WHITE_RESAMP_DSM.value % (self.trg_subject, atlas_suffix)))
        dipoles_files = [lh_dipoles_file, rh_dipoles_file]
        whites_resamp_dsm = [lh_white_dsm, rh_white_dsm]

        jobs1 = []
        jobs2 = []
        jobs3 = []
        jobs4 = [[], []]
        for ih, hemi in enumerate(["lh", "rh"]):
            jobs1.append(Job(ResamplingJobNames.MRI_SURF2SURF.value,
                             node_label="mri_surf2surf resample %s.white" % hemi))
            jobs1[-1].addArguments("--srcsubject", self.subject, "--trgsubject", self.trg_subject, "--hemi", hemi,
                              "--sval-xyz", "white", "--tval", "white-%s" % self.trg_subject, "--tval-xyz", t1_mgz)
            jobs1[-1].uses(t1_mgz, link=Link.INPUT)
            jobs1[-1].uses(whites[ih], link=Link.INPUT)
            jobs1[-1].uses(whites_resamp[ih], link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(jobs1[-1])

            dax.depends(jobs1[-1], job_head_model)

            jobs2.append(Job(HeadModelJobNames.CONVERT_TO_BRAIN_VISA.value,
                             node_label="%s resampled white surface conversion to brain visa" % hemi))
            jobs2[-1].addArguments(whites_resamp[ih], whites_resamp_tri[ih], self.subject)
            jobs2[-1].uses(whites_resamp[ih], link=Link.INPUT)
            jobs2[-1].uses(whites_resamp_tri[ih], link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(jobs2[-1])

            dax.depends(jobs2[-1], jobs1[-1])

            jobs3.append(Job(HeadModelJobNames.OM_ASSEMBLE.value,
                             node_label="om_assemble %s surface source model" % hemi))
            jobs3[-1].addArguments("-SurfSourceMat", head_model_geom, head_model_cond, whites_resamp_tri[ih],
                              whites_resamp_ssm[ih])
            for surf in bem_tri_surfs:
                jobs3[-1].uses(surf, link=Link.INPUT)
            jobs3[-1].uses(head_model_geom, link=Link.INPUT)
            jobs3[-1].uses(head_model_cond, link=Link.INPUT)
            jobs3[-1].uses(whites_resamp_tri[ih], link=Link.INPUT)
            jobs3[-1].uses(whites_resamp_ssm[ih], link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(jobs3[-1])

            dax.depends(jobs3[-1], jobs2[-1])

            for iatlas, atlas_suffix in enumerate(self.atlas_suffixes):

                if len(atlas_suffix) == 0:
                    atlas_name = "default"
                else:
                    atlas_name = atlas_suffix[1:]

                jobs4[ih].append(Job(HeadModelJobNames.OM_ASSEMBLE.value,
                                     node_label="om_assemble %s deep source model for atlas %s" % (hemi, atlas_name)))
                jobs4[ih][-1].addArguments("-DipSourceMat", head_model_geom, head_model_cond,
                                           dipoles_files[ih][iatlas], whites_resamp_dsm[ih][iatlas])
                for surf in bem_tri_surfs:
                    jobs4[ih][-1].uses(surf, link=Link.INPUT)
                jobs4[ih][-1].uses(head_model_geom, link=Link.INPUT)
                jobs4[ih][-1].uses(head_model_cond, link=Link.INPUT)
                jobs4[ih][-1].uses(dipoles_files[ih][iatlas], link=Link.INPUT)
                jobs4[ih][-1].uses(whites_resamp_dsm[ih][iatlas], link=Link.OUTPUT, transfer=True, register=True)
                dax.addJob(jobs4[ih][-1])

                dax.depends(jobs4[ih][-1], jobs_mapping_details[iatlas])
                dax.depends(jobs4[ih][-1], jobs3[-1])

        return jobs4
