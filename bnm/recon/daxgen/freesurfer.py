import os
import time
import sys
import enum

# XXX qualified name preferable
from Pegasus.DAX3 import *


class PropKeys(enum.Enum):
    "Enumerates possible keys in a properties file."
    subject = "subject"
    subject_folder = "subject.folder"
    threads = "open.mp.threads"
    t1_folder = "t1.folder"
    t1_input_format = "t1.format"
    t2_flag = "t2.flag"
    t2_folder = "t2.folder"
    mri_folder = "mri.folder"
    flair_flag = "flair.flag"
    flair_folder = "flair.folder"


class FileNames(enum.Enum):
    "Enumerates special FreeSurfer file names."
    mri_aparc_aseg = "aparc+aseg"


class FileTypes(enum.Enum):
    "Enumerates file types."
    dicom = "dicom"
    mgz = "mgz"
    nifti = "nifti"


# XXX make a cli arg
CONFIG_FILE = "config/patient.properties"


def read_from_properties_file():
    config = dict()
    with open(CONFIG_FILE) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or len(line) < 3 or line.index("=") < 0:
                continue
            key_str, value = line.split("=")
            try:
                key = PropKeys(key_str)
            except KeyError:
                raise Exception('Invalid property key %r in file %r.' % (
                        key_str, CONFIG_FILE))
            config[key] = value
    print "Read patient configuration", config
    return config


def generate_uq_file_name(parent_folder, name_pattern):
    file_name = name_pattern.replace("*", str(time.clock()))
    return os.path.join(parent_folder, file_name)


def step_recon_all_1(t1_in_path, t1_format, t1_out_path):
    # # T1 input pre-processing
    # if [ "$T1_INPUT_FRMT" = "dicom" ]
    # then
    #     mri_convert $T1 $T1/t1_raw.nii.gz --out_orientation RAS -rt nearest
    #     T1=$T1/t1_raw.nii.gz
    # fi
    # #ENDIF

    t1_in = File(t1_in_path)
    t1_out = File(t1_out_path)
    print "Processing for T1 format", t1_format
    if t1_format == FileTypes.dicom:
        print " - Dicom identified"
        job = Job(name="mri_convert", node_label="T1 input pre-processing")
        job.addArguments(t1_in, t1_out, "--out_orientation", "RAS", "-rt", "nearest")
    else:
        print " - Keep the original"
        job = Job(name="cp", node_label="Copy t1_in_path into t1_out_path")
        job.addArguments(t1_in, t1_out)

    job.uses(t1_in, link=Link.INPUT)
    job.uses(t1_out, link=Link.OUTPUT, transfer=False, register=False)
    return job


def step_recon_all_2(t1_in_path, subject_name, open_mp_threads):
    # # Freesurfer T1 processing
    # recon-all -s ${SUBJECT} -i $T1 -all -parallel -openmp $OPENMP_THRDS
    t1_in = File(t1_in_path)
    job = Job(name="recon-all", node_label="Call 'recon-all'")
    job.addArguments("-s", subject_name, "-i", t1_in, "-all", "-parallel", "-openmp", open_mp_threads)
    job.uses(t1_in, link=Link.INPUT)
    return job


def step_recon_all_3(mri_folder, mri_out_path):
    # # This could be made already here:
    # # Generate nifti file with good orientation
    # mri_convert $MRI/aparc+aseg.mgz $MRI/aparc+aseg.nii.gz --out_orientation RAS -rt nearest
    mri_in = File(os.path.join(mri_folder, FileNames.mri_aparc_aseg.value + ".mgz"))
    mri_out = File(mri_out_path)
    job = Job(name="mri_convert", node_label="Generate nifti file with good orientation")
    job.addArguments(mri_in, mri_out, "--out_orientation", "RAS", "-rt", "nearest")
    job.uses(mri_in, link=Link.INPUT)
    job.uses(mri_out, link=Link.OUTPUT, transfer=False, register=False)
    return job


def step_snapshot():
    # # Visual checks (screen-shots) for brain, pial, white,
    # bash snapshot.sh white+pial $SUBJ_DIR/mri/T1.mgz $SUBJ_DIR/surf
    # for h in lh rh
    # do
    # 	bash snapshot.sh surf $SUBJ_DIR/surf/$h.inflated aparc
    # done
    job = Job(name="snapshot", node_label="Generate visual checks for this step")
    # TODO
    return job


def freesurfer_main():
    # The name of the DAX file is the first argument
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: %s DAXFILE\n" % (sys.argv[0]))
        sys.exit(1)
    daxfile = sys.argv[1]

    dax = ADAG("BrainNetworkModelReconstructionFlow")

    if hasattr(dax, 'metadata'):
        dax.metadata("created", time.ctime())
    config = read_from_properties_file()

    path_t1_raw = config[PropKeys.t1_folder]
    path_t1 = generate_uq_file_name(path_t1_raw, "t1_*.nii.gz")

    step1 = step_recon_all_1(path_t1_raw, config[PropKeys.t1_input_format], path_t1)
    dax.addJob(step1)

    step2 = step_recon_all_2(path_t1, config[PropKeys.subject], config[PropKeys.threads])
    dax.addJob(step2)

    path_mri_folder = config[PropKeys.mri_folder]
    path_mri = generate_uq_file_name(path_mri_folder, FileNames.mri_aparc_aseg.value + "_*.nii.gz")
    step3 = step_recon_all_3(path_mri_folder, path_mri)
    dax.addJob(step3)

    # Add control-flow dependencies
    dax.depends(step2, step1)
    dax.depends(step3, step2)

    # Write the DAX to stdout
    print "Writing %s" % daxfile
    f = open(daxfile, "w")
    dax.writeXML(f)
    f.close()
