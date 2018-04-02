# Reconstruction pipeline

This pipeline is offering a solution to build full brain network models starting from standard structural MR scans.
It is used to preprocess the MR scans in order to get actual files that are compatible with TVB. The result can be later uploaded in TVB. 

The mandatory inputs are DWI and T1 scans. Optionally, CT scans can be given as input, if sensors preprocessing is needed.

We are using the Pegasus WMS in order to connect and automatize the pipeline steps.

## Folder structure
- data
	<br> Here is where we keep some example data files. Some of them are intermediate files generated during pipeline run. There is also a minimal set of files that defines a TVB head. These files are currently used only for tests.  
- docs
	<br> This folder holds some visual hints of the dependencies between the pipeline steps. For example, there is an overview diagram representing the pipeline stages. These stages are defined by scripts that can be found inside the bin folder. On the other hand, there is also an example graph diagram which is displaying the more detailed steps. This kind of diagram is automatically generated at each pipeline run.  
- pegasus
	<br> All the configuration files necessary for the software run are kept here. 
	<br> At the first level there are the entry points which are explained bellow in the <b><i>Entry Point</i></b> section. 
	<br> Inside the <i>config</i> folder, there is a folder called <i>scripts</i>. This holds bash files with calls to pipeline commands that need some specific environment configuration. These bash files should not be changes, neither their names, because they are mapped inside the workflow.
	<br> The config folder also contains the configuration files which are patient specific. These are the only files that will change for each run and/or patient. The possible configurations are explained in a later section.
- provision
	<br> Inside this folder, you can find more details about the project dependencies.
- tvb
	<br> Actual python implementation is kept here. 
	<br> Package <i>dax</i> is using the Pegasus API to define a pipeline workflow and generate the jobs graph.
	<br> Inside <i>algo</i> package we have all the computations, services and algorithms used during the pipeline steps. Here we also have the <i>reconutils.py</i> which defines an API for calls made from the <i>config/scripts</i> bash files.
	<br> The model classes can be found inside <i>model</i> package.
	<br> Package <i>io</i> provides read/write functionalities from/to a variety of file formats.
	<br> We also have a set of tests for this module and they are kept inside the <i>tests</i> package.

## How to launch

### Entry point:
There are 2 available entry points for the pipeline. They are both under <i>pegasus</i> folder. In order to use these entry points, there are, in both cases, some configurations to be defined first. These configurations are kept as a folder specific to each patient and are explained in the next section.

The pipeline can be started using one of the following entry points:
- <b>main_pegasus.sh</b>
    <br> This is the most straight-forward one. It starts one pipeline run for a single patient based on a set of predefined configurations.
    <br> Command to launch the pipeline with this script: <b><i>sh main_pegasus.sh path_to_configurations_folder path_to_dax_folder</i></b>
    <br> The argument <i>path_to_configurations_folder</i> represents the path to the patient configuration files and <i>path_to_dax_folder</i> represents the folder where the dax will be generated.

- <b>run_sequential.py</b> 
    <br> This is a little more complex. It is used to start pipeline runs for a list of patients with similar configurations. As the name is suggesting, the runs will be started sequentially for the patients.
    <br> Command to launch the pipeline with this script: <b><i>python run_sequentially.py</i></b>
    <br> This script does not need arguments, but it needs the user to define the necessary configurations inside it. The necessary configurations are described bellow:
    - <b>PATH_TO_INPUT_SUBJ_FOLDERS</b>: path to the folder where you keep your patient raw data 
    - <b>PATH_TO_SUBJ_CONFIG_FOLDERS</b>: path to the folder where you keep your patient configurations
    - <b>PATH_TO_OUTPUT_SUBJ_FOLDER</b>: path to the folder where you want your outputs to be saved
    - <b>PREFIX_SUBJECT_FOLDER</b>: prefix of the patient folder name
    - <b>SUBJECTS_TO_BE_PROCESSED</b>: suffix of the patient folder name
    - <b>ATLASES</b>: a list with all the atlases to be used
    - <b>PATH_TO_DEFAULT_PEGASUS_CONFIGURATION</b>: path to the folder where you keep the default configuration files. E.g: tvb-recon/pegasus/config
    
### Configurations
All the configuration files are under <i>pegasus/config</i> at the top level. There are configurations specific to the patient, to the machine where the workflow is running or to the actual run. Some details about each file, are given bellow:
- environment_config.sh
    <br> These are the machine specific configurations. The configurations are mostly paths of the software and other variables that are needed for the pipeline environment setup.
- patient_flow.properties
    <br> These are the patient specific configurations. 
- pegasus.properties
    <br> This is a pegasus specific configuration file. It just defined the paths to sites.xml, rc.txt, tc.txt and rc_out.txt.
- rc.txt
    <br> Inside this file, the inputs should be defined in a key-value format, where value is the path to the input file.
- rc_out.txt
    <br> Using this file, the output can be structured in a similar way as the input.
- sites.xml
    <br> This should not change since the configurable variables will be taken from the environment.
- tc.xml
    <br> This should not be change since it contains the commands mapping.


## How to rerun

### Rerun with different parameters

### Recover over shutdown


## Environment
The pipeline steps are dependent on the following external tools:
- <b>Freesurfer</b>
	<br> Freesurfer is used to process the T1 images and generate the anatomy of a patient brain.
- <b>Mrtrix</b>
	<br> Mrtrix is used to preprocess the DWI images and also to generate the tracts.
- <b>FSL</b>
	<br> FSL is used to register images from different spaces to a common one.
- <b>MNE</b>
	<br> MNE is used to generate BEM surfaces.
- <b>OpenMEEG</b>
	<br> OpenMEEG can be used for sensor computations.

The automatized workflow is based on:
- <b>Pegasus</b>
	<br> This is the workflow engine we have used for automatizing the pipeline steps.
- <b>HTCondor</b>
	<br> Pegasus uses HTCondor as a job scheduler.

[![Build Status](https://travis-ci.org/the-virtual-brain/tvb-recon.svg?branch=master)](https://travis-ci.org/the-virtual-brain/tvb-recon) [![Coverage Status](https://coveralls.io/repos/github/the-virtual-brain/tvb-recon/badge.svg)](https://coveralls.io/github/the-virtual-brain/tvb-recon)
