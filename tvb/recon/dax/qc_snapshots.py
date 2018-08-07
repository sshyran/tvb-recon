from Pegasus.DAX3 import ADAG, Job, Link, File


class QCSnapshots(object):
    __instance = None
    SNAPSHOT_NUMBER = 0

    def __init__(self):
        if self.__instance is not None:
            raise ValueError("This is a singleton and an instance already exists!")

    @classmethod
    def get_instance(cls):
        if cls.__instance is None:
            cls.__instance = QCSnapshots()
        return cls.__instance

    def construct_basename(self, name=""):
        basename = "snapshot"
        if len(name) > 0:
            basename = basename + "_" + name
        return basename

    def add_2vols_snapshot_step(self, dax, jobs_before, vol1, vol2, name=""):
        basename = self.construct_basename(name)
        snapshot_file_1 = File(basename + "_sagittal_%d.png" % self.SNAPSHOT_NUMBER)
        snapshot_file_2 = File(basename + "_coronal_%d.png" % self.SNAPSHOT_NUMBER)
        snapshot_file_3 = File(basename + "_axial_%d.png" % self.SNAPSHOT_NUMBER)

        job = Job("qc_snapshot", node_label=basename)
        job.addArguments(str(self.SNAPSHOT_NUMBER), basename, "2vols", vol1, vol2)
        job.uses(vol1, link=Link.INPUT)
        job.uses(vol2, link=Link.INPUT)

        job.uses(snapshot_file_1, link=Link.OUTPUT, transfer=True, register=False)
        job.uses(snapshot_file_2, link=Link.OUTPUT, transfer=True, register=False)
        job.uses(snapshot_file_3, link=Link.OUTPUT, transfer=True, register=False)
        dax.addJob(job)

        for job_before in jobs_before:
            dax.depends(job, job_before)

        self.SNAPSHOT_NUMBER += 1

    def add_3vols_snapshot_step(self, dax, jobs_before, vol1, vol2, vol3, name=""):
        basename = self.construct_basename(name)
        snapshot_file_1 = File(basename + "_sagittal_%d.png" % self.SNAPSHOT_NUMBER)
        snapshot_file_2 = File(basename + "_coronal_%d.png" % self.SNAPSHOT_NUMBER)
        snapshot_file_3 = File(basename + "_axial_%d.png" % self.SNAPSHOT_NUMBER)

        job = Job("qc_snapshot", node_label=basename)
        job.addArguments(str(self.SNAPSHOT_NUMBER), basename, "3vols", vol1, vol2, vol3)
        job.uses(vol1, link=Link.INPUT)
        job.uses(vol2, link=Link.INPUT)
        job.uses(vol3, link=Link.INPUT)

        job.uses(snapshot_file_1, link=Link.OUTPUT, transfer=True, register=False)
        job.uses(snapshot_file_2, link=Link.OUTPUT, transfer=True, register=False)
        job.uses(snapshot_file_3, link=Link.OUTPUT, transfer=True, register=False)
        dax.addJob(job)

        for job_before in jobs_before:
            dax.depends(job, job_before)

        self.SNAPSHOT_NUMBER += 1

    def add_vol_surf_snapshot_step(self, dax, jobs_before, vol, surfs, name=""):
        basename = self.construct_basename(name)
        snapshot_file_1 = File(basename + "_sagittal_%d.png" % self.SNAPSHOT_NUMBER)
        snapshot_file_2 = File(basename + "_coronal_%d.png" % self.SNAPSHOT_NUMBER)
        snapshot_file_3 = File(basename + "_axial_%d.png" % self.SNAPSHOT_NUMBER)

        job = Job("qc_snapshot", node_label=basename)
        job.addArguments(str(self.SNAPSHOT_NUMBER), basename, "vol_surf", vol)
        for surf in surfs:
            job.addArguments(surf)
        job.uses(vol, link=Link.INPUT)
        for surf in surfs:
            job.uses(surf, link=Link.INPUT)
        job.uses(snapshot_file_1, link=Link.OUTPUT, transfer=True, register=False)
        job.uses(snapshot_file_2, link=Link.OUTPUT, transfer=True, register=False)
        job.uses(snapshot_file_3, link=Link.OUTPUT, transfer=True, register=False)
        dax.addJob(job)

        for job_before in jobs_before:
            dax.depends(job, job_before)

        self.SNAPSHOT_NUMBER += 1

    def add_surf_annot_snapshot_step(self, dax, jobs_before, surf, annot, name=""):
        basename = self.construct_basename(name)
        snapshot_file_1 = File(basename + "_surface_annotation_%d0.png" % self.SNAPSHOT_NUMBER)
        snapshot_file_2 = File(basename + "_surface_annotation_%d1.png" % self.SNAPSHOT_NUMBER)
        snapshot_file_3 = File(basename + "_surface_annotation_%d2.png" % self.SNAPSHOT_NUMBER)
        snapshot_file_4 = File(basename + "_surface_annotation_%d3.png" % self.SNAPSHOT_NUMBER)
        snapshot_file_5 = File(basename + "_surface_annotation_%d4.png" % self.SNAPSHOT_NUMBER)
        snapshot_file_6 = File(basename + "_surface_annotation_%d5.png" % self.SNAPSHOT_NUMBER)

        job = Job("qc_snapshot", node_label=basename)
        job.addArguments(str(self.SNAPSHOT_NUMBER), basename, "surf_annot", surf, annot)
        job.uses(surf, link=Link.INPUT)
        job.uses(annot, link=Link.INPUT)

        job.uses(snapshot_file_1, link=Link.OUTPUT, transfer=True, register=True)
        job.uses(snapshot_file_2, link=Link.OUTPUT, transfer=True, register=True)
        job.uses(snapshot_file_3, link=Link.OUTPUT, transfer=True, register=True)
        job.uses(snapshot_file_4, link=Link.OUTPUT, transfer=True, register=True)
        job.uses(snapshot_file_5, link=Link.OUTPUT, transfer=True, register=True)
        job.uses(snapshot_file_6, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job)

        for job_before in jobs_before:
            dax.depends(job, job_before)

        self.SNAPSHOT_NUMBER += 1
