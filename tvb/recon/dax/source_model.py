from Pegasus.DAX3 import Job, File, Link

from tvb.recon.dax.mappings import T1Files, ResamplingJobNames, SourceModelFiles, HeadModelJobNames, AsegFiles, \
    HeadModelFiles


class SourceModel(object):
    def __init__(self, subject, trg_subject):
        self.subject = subject
        self.trg_subject = trg_subject

    def add_source_model_steps(self, dax, job_head_model, job_mapping_details):
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

        lh_dipoles_file = File(AsegFiles.LH_DIPOLES_TXT.value)
        rh_dipoles_file = File(AsegFiles.RH_DIPOLES_TXT.value)
        dipoles_files = [lh_dipoles_file, rh_dipoles_file]

        lh_white_dsm = File(SourceModelFiles.LH_WHITE_RESAMP_DSM.value % self.trg_subject)
        rh_white_dsm = File(SourceModelFiles.RH_WHITE_RESAMP_DSM.value % self.trg_subject)
        whites_resamp_dsm = [lh_white_dsm, rh_white_dsm]

        last_job = None

        for idx, hemi in enumerate(["lh", "rh"]):
            job1 = Job(ResamplingJobNames.MRI_SURF2SURF.value)
            job1.addArguments("--srcsubject", self.subject, "--trgsubject", self.trg_subject, "--hemi", hemi,
                              "--sval-xyz", "white", "--tval", "white-%s" % self.trg_subject, "--tval-xyz", t1_mgz)
            job1.uses(t1_mgz, link=Link.INPUT)
            job1.uses(whites[idx], link=Link.INPUT)
            job1.uses(whites_resamp[idx], link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(job1)

            dax.depends(job1, job_head_model)

            job2 = Job(HeadModelJobNames.CONVERT_TO_BRAIN_VISA.value)
            job2.addArguments(whites_resamp[idx], whites_resamp_tri[idx], self.subject)
            job2.uses(whites_resamp[idx], link=Link.INPUT)
            job2.uses(whites_resamp_tri[idx], link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(job2)

            dax.depends(job2, job1)

            job3 = Job(HeadModelJobNames.OM_ASSEMBLE.value)
            job3.addArguments("-SurfSourceMat", head_model_geom, head_model_cond, whites_resamp_tri[idx],
                              whites_resamp_ssm[idx])
            for surf in bem_tri_surfs:
                job3.uses(surf, link=Link.INPUT)
            job3.uses(head_model_geom, link=Link.INPUT)
            job3.uses(head_model_cond, link=Link.INPUT)
            job3.uses(whites_resamp_tri[idx], link=Link.INPUT)
            job3.uses(whites_resamp_ssm[idx], link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(job3)

            dax.depends(job3, job2)

            job4 = Job(HeadModelJobNames.OM_ASSEMBLE.value)
            job4.addArguments("-DipSourceMat", head_model_geom, head_model_cond, dipoles_files[idx],
                              whites_resamp_dsm[idx])
            for surf in bem_tri_surfs:
                job4.uses(surf, link=Link.INPUT)
            job4.uses(head_model_geom, link=Link.INPUT)
            job4.uses(head_model_cond, link=Link.INPUT)
            job4.uses(dipoles_files[idx], link=Link.INPUT)
            job4.uses(whites_resamp_dsm[idx], link=Link.OUTPUT, transfer=True, register=True)
            dax.addJob(job4)

            dax.depends(job4, job_mapping_details)
            dax.depends(job4, job3)

            last_job = job4

        return last_job
