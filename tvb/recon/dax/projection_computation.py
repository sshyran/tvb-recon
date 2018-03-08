from Pegasus.DAX3 import File, Job, Link
from tvb.recon.dax.mappings import ProjectionCompFiles, ProjectionCompJobNames, AsegFiles


class ProjectionComputation(object):
    def __init__(self, subject, sensors_type, atlas_suffix):
        self.subject = subject
        self.sensors_type = sensors_type
        self.atlas_suffix = atlas_suffix

    def add_projection_computation_steps(self, dax, job_mapping_details):
        projection_mat = File(ProjectionCompFiles.PROJECTION_MAT.value % (self.sensors_type, self.atlas_suffix))
        sensor_positions = File(ProjectionCompFiles.SENS_POSITIONS.value % self.sensors_type)
        centers_txt = File(AsegFiles.CENTERS_TXT.value % self.atlas_suffix)
        job = Job(ProjectionCompJobNames.COMPUTE_PROJ_MAT.value)
        job.addArguments(sensor_positions, centers_txt, projection_mat, self.subject)
        job.uses(sensor_positions, link=Link.INPUT)
        job.uses(centers_txt, link=Link.INPUT)
        job.uses(projection_mat, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job)

        dax.depends(job, job_mapping_details)
