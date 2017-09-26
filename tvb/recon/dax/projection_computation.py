from Pegasus.DAX3 import File, Job, Link
from tvb.recon.dax.mappings import ProjectionCompFiles, ProjectionCompJobNames, OutputConvFiles


class ProjectionComputation(object):
    def __init__(self, subject, sensors_type):
        self.subject = subject
        self.sensors_type = sensors_type

    def add_projection_computation_steps(self, dax, job_conn):
        projection_mat = File(ProjectionCompFiles.PROJECTION_MAT.value % self.sensors_type)
        sensor_positions = File(ProjectionCompFiles.SENS_POSITIONS.value % self.sensors_type)
        connectivity_zip = File(OutputConvFiles.CONNECTIVITY_ZIP.value)
        job = Job(ProjectionCompJobNames.COMPUTE_PROJ_MAT.value)
        job.addArguments(sensor_positions, connectivity_zip, projection_mat, self.subject)
        job.uses(sensor_positions, link=Link.INPUT)
        job.uses(connectivity_zip, link=Link.INPUT)
        job.uses(projection_mat, link=Link.OUTPUT, transfer=True, register=True)
        dax.addJob(job)

        dax.depends(job, job_conn)
