# Reconstruction tools

Tools for building full brain network models from standard structural MR scans.

## Goals

- Scripts to help set up an environment with dependecies
- Extensive visual diagnostics for intermediate & final results
- Functionality for rapidly reparcellating and obtaining new brain model
- Adapters & datatypes for TVB framework

## TODO

- all subjects in same space & orientation
- use a single template electrode setup
- consider 128 or 256 EEG sensor set
- MEG should just work once template adjusted
- EEG shoudl move sensors to be a few mm outside the head surface (lh.seghead)
- use pg tool to adjust template layout against big head

- subparcellations via ?h.sphere (not volumetric)
- compute richer stats on connectome (tract len dist, FA dist)

FS has many tricks. It is worth tab-completing `mri_` and `mris_` to see what
stuff is possible.

- `mris_convert --combinesurfs <in1> <in2> <out>` combines surface files
- `mris_anatomical_stats -a label/lh.aparc.annot -b tvb lh` per label stats


## CGAL

Considered for

- [improved mesh simplification](http://doc.cgal.org/latest/Surface_mesh_simplification/)
- [geodesic distance algorithms](http://doc.cgal.org/latest/Surface_mesh_shortest_path/index.html)

```
sudo dnf install CGAL-devel CGAL
```

One potential advantage of CGAL's mesh simplification is that arbitrary user data can
be attached and used during simplification to maintain for example the surfacic
parcellation.

## outline steps

- fs recon all
- resample surface
- coreg t1 dwi
- dmr estimates
- tractography
- em fwd for cortical surface
- coreg t1 ct
- label electrodes

then for every connectivity, generate

- connectome matrices
- region mapping
- reduced lead fields

try to follow FS' organization in per-subject topic folders

## parameters

Most steps are parameter free, but a few to keep in mind

- csd max harmonic, lmax (can just be left at default, max order possible given data)
- subparcellation
- number of tracks (though 10M sifted to 5M seems safe default)

Act & sift always on, topup depends on data.

## new parcellations

- new subdivided parcellations can be made on a sphere
- save to new annot files label/{r,l}h.foo.annot
- mri_aparc2aseg --s tvb --aseg aseg.presurf.hypos --annot foo
- results in mri/foo+aseg.mgz
- usuable for tck2connectome

it seems nibabel is not writing good annotations, we'll need a fix
or workaround..


## From Scratch

- setting up a VM with Cent OS 7
- scripts to handle openmeeg, mrtrix compilation
- each step with example results & time estimates

- freesurfer -> bem/headmat, dmr/trac, decim surf, ct/elec
- subparc & trac -> conn
- conn & dceim surf & headmat -> eeg/meg/bold fwd
- & ct/elec -> seeg fwd

aparc+aseg is a spatial filter for generated
BOLD and should be considered essential for that.

### reproducible workflow

It's always a challenge to make a complex software workflow reproducible. In the
following, we're going to assume a CentOS 7 64-bit virtual machine, based off the
minimal installer. It's hoped that these instructions will allow anyone regardless
of platform (Windows, Mac, etc) to take advantage of the virtualization workflow
in a standardized fashion.

Various times are given below for running the workflow for a single subject at
a time, and for reference, the host CPU is a Core i7
2620M, the host HD is an Intel 530 Series SSD. The guest VM is CentOS 7 64-bit,
running in QEMU with KVM acceleration with 1 GB RAM & 30 GB disk. Other virtualization
software such as VMWare, Xen, VirtualBox, Parallels, etc should work fine too.

Because the setup procedure for a VM varies between software, we won't provide
instructions here for this, as they widely available on the internet.

### installing software

#### Setting up the system

Assuming your CentOS 7 64-bit VM is freshly set up, you should log
in as root, and download the scripts accompanying this guide, and run
the system setup script, which will update the package list, install required
pacakges and reboot
```bash
yum install -y git
export TVBVIRT=$HOME/tvbvirt
git clone https://github.com/maedoc/tvb-virtualizer $TVBVIRT
source $TVBVIRT/setup-system-and-reboot.sh
```
This will reboot the VM into a graphical environment. Log in as root (click "Not Listed?",
and type `root` as username), then open a Terminal (Applications > Utilities > Terminal)
in order to continue.

The font in the Terminal looks oddly spaced, to fix, in Edit > Profile Preferences,
check Custom font and select Liberation Mono Regular.

#### Build MRtrix3 & OpenMEEG from source

These two tools must be built from source, so in the terminal, run the provided
build script:
```bash
export TVBVIRT=$(pwd)/tvbvirt
source $TVBVIRT/setup-software.sh
```
This takes about X minutes.

Finally, depending on your system's configuration, you may wish to limit
the number of threads used by MRtrix. By default, only one core is used,
though you may wish to edit `/etc/mrtrix.conf` to set this and
other options as described on the
[MRtrix configuration file wiki page](https://github.com/MRtrix3/mrtrix3/wiki/MRtrix-configuration-file).

#### OpenGL 3.3 support for MRtrix GUIs

Compile a newer version of Mesa3d & llvmpipe driver,
_TODO_.

#### FreeSurfer

Setting up FreeSurfer requires
- Downloading it
- Unpacking it
- Registering
- Placing your license file in the installation

The first two steps are handled by a provided script,
```bash
source $TVBVIRT/setup/50-freesurfer.sh
```
This can take a while as FreeSurfer is a large download.

While it's downloading, you should [complete registration](https://surfer.nmr.mgh.harvard.edu/fswiki/Registration)
(which is free), and they will send you a small license file with
some text like
```
foo@bar.com
123456
 *jf0a0fj2932j
 jsrl900r9u309
```
which you should place in the license file in the FreeSurfer folder,
```bash
cat > /usr/local/freesurfer/license.txt <<EOF
foo@bar.com
123456
 *jf0a0fj2932j
 jsrl900r9u309
EOF
```

#### FSL

FSL provides Python script for installation on Linux, and you are required to
- open Firefox in the VM (Applications > Internet > Firefox)
- go to their site to agree to their license: http://fsl.fmrib.ox.ac.uk/fsldownloads
- enter your information
- in section `4. SELECT AND DOWNLOAD`, click the `fslinstaller.py` button.

This results in a file `Downloads/fslinstaller.py`. In the Terminal,
```bash
python Downloads/fslinstaller.py
```
Accept the default installation location (`/usr/local`), and wait for it to finish downloading.

If the download fails, check the size of your tmpfs; if it is too small for the FSL
download, this can cause the download to fail, and you can change in the `fslinstaller.py`
scripts `mkstemp()` call to `mkstemp(dir=os.getcwd())`, and re-run the installer.

#### MNE tools

MNE tools are yet another package which requires manual registration and download from

http://www.nmr.mgh.harvard.edu/martinos/userInfo/data/MNE_register/index.php

After registration, enter your email and password, and click `Download`. On the
resulting page, click the `MNE-2.7.0-3106-Linux-x86_64.tar.gz` link.

Once the download is finished, in the Terminal,

```bash
source $TVBVIRT/setup/60-mne.sh
```

#### Setting up a default environment

By default FreeSurfer saves subject data inside its own subjects directory
but this is not writeable except by root, and in general it's a good
idea to specify your own subjects dir, so in your regular user (not root,
so after all the above set up, log out and log back in as a non-root user)
home directory,
```
export TVBVIRT=$HOME/tvbvirt
git clone https://github.com/maedoc/tvb-virtualizer $TVBVIRT
source $TVBVIRT/setup-user-env.sh
```
This will create a subjects folder and add environment variables to allow
all the different software installed to function correctly.

Close the Terminal and open a new one. At this point, all commands used in the
virtualization workflow should function.

## Workflow

In the following sections we step through the workflow used for the
default dataset. Note that the dMR data were already corrected for
eddy distortion with FSL's `topup`, but a script is provided to handle
that if you need it.

This workflow produces

1. Cortical surface usuable for simulation
2. ROI diffusion based connectome, usuable for simulation
3. Cortical parcellation for mapping cortical regions in (2) to vertices in (1)
4. Volumetric parcellation for cortical & subcortical regions in (2)
5. BEM surfaces, including a scalp/head surface
6. Forward solutions from activty on ROIs in (2) to sEEG, EEG & MEG, w/ sensors xyz & ori
7. BOLD forward solution via (2)
8. Cortical surface geodesic distance matrix

The dependencies are, roughly,

1. T1
2. DMR
3. CT
4. recon-all -> 1
5. Low-res surface & annot -> 4


Producing these data requires many intermediate data that will be organized
in the subject's folder, following FreeSurfer conventions as possible.

Finally, in the following examples, the subject name is always `tvb`, so
by convention we set this in an environment variable before starting:
```bash
export SUBJECT=tvb
```

### FreeSurfer recon-all

```bash
recon-all -s ${SUBJECT} -i t1_raw.nii.gz -all
```
This takes takes 11 hours, though the time required can vary
significantly across datasets.

Loading the T1, aparc+aseg, lh.white surface with lh.aparc.DK annot should look
like

![this](img/recon-all.png)

### Higher resolution parcellations

The default parcellation can be subdivided into regions of roughly 250.0 mm^2
```bash
mris_divide_parcellation ${SUBJECT} lh aparc 250.0 aparc250
mris_divide_parcellation ${SUBJECT} rh aparc 250.0 aparc250
```
This takes about 10 seconds. Check the subdivision visually,
```
freeview -f ${SUBJECTS_DIR}/${SUBJECT}/surf/lh.pial:annot=aparc250 -viewport 3d
```
![aparc250](img/aparc250.png)

For each such parcellation, the corresponding volume should be generated
```
mri_aparc2aseg --s ${SUBJECT} --aseg aseg --annot aparc250
```
This takes 130 seconds and generates `mri/aparc250+aseg.mgz`.
```
freeview -v ${SUBJECTS_DIR}/${SUBJECT}/mri/aparc250+aseg.mgz -viewport 3d
```
![aseg250](img/aparc250-aseg.png)

### Lower resolutions surfaces

The cortical surfaces generated by FreeSurfer are around 150k vertices, which
is far too many for simulation on current hardware, so we'll map the surface and
its parcellation to a lower resolution (~10k) using the `fsaverage5` subjects,

```bash
cp -r ${FREESURFER_HOME}/subjects/fsaverage5 ${SUBJECTS_DIR}
mri_surf2surf --srcsubject ${SUBJECT} --trgsubject fsaverage5 --hemi lh --sval-xyz pial --tval pial.${SUBJECT} --tval-xyz ${SUBJECTS_DIR}/${SUBJECT}/mri/T1.mgz
mri_surf2surf --srcsubject ${SUBJECT} --trgsubject fsaverage5 --hemi rh --sval-xyz pial --tval pial.${SUBJECT} --tval-xyz ${SUBJECTS_DIR}/${SUBJECT}/mri/T1.mgz
cp ${SUBJECTS_DIR}/fsaverage5/surf/lh.pial.${SUBJECT} ${SUBJECTS_DIR}/${SUBJECT}/surf/lh.pial.low
cp ${SUBJECTS_DIR}/fsaverage5/surf/rh.pial.${SUBJECT} ${SUBJECTS_DIR}/${SUBJECT}/surf/rh.pial.low
```
produces ~10k vertex surfaces. Inspect the result visually:
```bash
freeview -v ${SUBJECTS_DIR}/${SUBJECT}/mri/T1.mgz -f ${SUBJECTS_DIR}/${SUBJECT}/surf/*.pial.low -viewport coronal
```
![this](img/t1-pial-low.png)

_This seems to work only with the dev version of FS, looking for a workaround
to embed source geometry in target surface, on FS stable._

Next, we obtain the parcellations for the lower resolution surfaces,

```bash
for parc in aparc aparc250
do
  for h in lh rh
  do
    mri_surf2surf --srcsubject ${SUBJECT} --trgsubject fsaverage5 --hemi $h \
      --sval-annot ${SUBJECTS_DIR}/${SUBJECT}/label/$h.$parc.annot \
      --tval ${SUBJECTS_DIR}/${SUBJECT}/label/$h.low.$parc
  done
done
```
This takes about 10 seconds per parcellation, and the `aparc` should look like
![this](img/t1-pial-low-annot.png)

### Surface geodesic distance matrix

The geodesic matrices of vertex-vertex distances along the cortical
surfaces are computed as
```bash
python compute_gdist.py
```
This requires 150 seconds.

### sEEG sensor identification

Coregister the CT and T1 image
```bash
regopt="-dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -cost mutualinfo"
flirt -in ct.nii -ref t1.nii -omat ct2t1.mat -out ct-in-t1.nii.gz $regopt
```
This takes about 9 minutes for this dataset.
Checking the results visually
```bash
freeview -v t1.nii.gz ct.nii.gz:colormap=jet
```
![t1-ct](img/t1-ct.png)

The electrodes must be manually labeled, so keep Freeview open, open
whatever implantation reference availabe (x-ray, schema, etc), and
for each electrode, create a control point set with the name of the electrode
with two points 1) at the deepest part of the electrode and 2) where the
electrode intersects the skull. Make sure to save these point sets in the
`${SUBJECTS_DIR}/${SUBJECT}/sensors/seeg` directory.

Once all electrodes are labeled, generate the list of contacts
```bash
python gen_seeg_contacts.py
```
If you have the list of x-y-z positions of contacts by some other means,
then you can skip manual labelling and place this list in the
`${SUBJECTS_DIR}/${SUBJECT}/sensors/seeg-contacts.txt` file.
_Format to be defined._

### s/M/EEG forward solutions

#### BEM Surfaces

```
mne_watershed_bem --subject tvb
```
This takes about 12 minutes, and produces four surface files
`subjects/tvb/bem/watershed/tvb_{brain,inner_skull,outer_skin,outer_skull}_surface`.

Each has 10242 vertices, which is too high for further processing, so we'll
decimate them:

```bash
for surf in subjects/tvb/bem/watershed/*_surface
do
    time mris_decimate -d 0.1 ${surf} ${surf}-low
done
```
This takes about 10 seconds for the 4 surfaces.
Inspect the resulting BEM surfaces on the T1:
```
freeview -v subjects/tvb/mri/T1.mgz -f subjects/tvb/bem/watershed/*-low -viewport coronal
```
![this](img/t1-bem.png)

If for memory reasons you wish to decimate the surfaces further, the `mrsi_decimate_gui`
command will allow you to interactively decimate and visualize the results.

#### Head model

OpenMEEG doesn't currently read FreeSurfer format surfaces as produced by the
preview BEM Surfaces step, so we first convert them to BrainVisa format:
```
python convert_bem_to_tri.py
```
Now we can use the generic head model geometry and conductivity files
to generate & invert the head model:
```bash
python gen_head_model.py
pushd ${SUBJECTS_DIR}/${SUBJECT}/bem
om_assemble -HM head_model.geom head_model.cond head.mat # 2m32s
om_minverser head.mat head-inv.mat # 3m30s
popd
```
This requires about 6 minutes.

#### Gain matrices

- combine lh & rh surfaces
- for each sensor set
  - compute full surface fwd
  - for each parc
    - average or svd the ROI fwds

Subcortical structures however are better modeled as unconstrained
volumetric source space, though we should assess point spread & cross
talk to set grid size; so we need some work here, cf.

- http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0059856

Modeling volumetrically would mean

- every grid point has three dipoles, three neural masses
- they can be coupling by mean field or ?

### dMR

Preparing a diffusion MR sequence for tractography requires a few preprocessing
step. We'll use a new subject subfolder, and change into that folder for this
section of the workflow
```bash
mkdir -p ${SUBECTS_DIR}/${SUBJECT}/dmr
pushd ${SUBJECTS_DIR}/${SUBJECT}/dmr
```

#### Conversion

Copy the diffusion into this folder as `dwi.nii.gz`, along with the `bvecs` and `bvals`,
and convert the data to MRtrix format, extracting the brain mask and low-b image.
```bash
mrconvert dwi.nii.gz dwi.mif -fslgrad bvecs bvals
dwi2mask dwi.mif mask.mif
dwiextract dwi.mif lowb.mif -bzero
mrconvert lowb.mif lowb.nii.gz
```
These steps take about 90 seconds depending on file system and hard disk speed,
and can require several GB of RAM depending on the dataset.

#### Registration with T1

In order to apply the anatomical labeling obtained in previous steps to the
diffusion images and tractography, we need to convert the relevant images
and obtain a transform from the diffusion image coordinate system to that of the T1:

```bash
for image in T1 aparc+aseg aparc250+aseg
do
  mri_convert ../mri/$image.mgz $image.nii.gz --out_orientation RAS -rt nearest
done
for parc in aparc aparc250
do
  fslreorient2std $parc+aseg.nii.gz $parc+aseg-reo.nii.gz
  mv $parc+aseg-reo.nii.gz $parc+aseg.nii.gz
done
regopt="-dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -cost mutualinfo"
flirt -in lowb.nii.gz -ref T1.nii.gz -omat d2t.mat -out lowb-in-t1.nii.gz $regopt
```
This takes about 90 seconds. Next, we apply the inverse transform to the T1 and
the volumetric parcellation:
```bash
convert_xfm -omat t2d.mat -inverse d2t.mat
for parc in aparc aparc250
do
  flirt -applyxfm -in $parc+aseg.nii.gz -ref lowb.nii.gz \
        -out $parc+aseg-in-d.nii.gz -init t2d.mat -interp nearestneighbour
done
flirt -applyxfm -in T1.nii.gz -ref lowb.nii.gz \
    -out t1-in-d.nii.gz -init t2d.mat -interp nearestneighbour
```
This takes a few seconds. Visualize the result:
```bash
freeview -v t1-in-d.nii.gz lowb.nii.gz:colormap=heat aparc+aseg-in-d.nii.gz:colormap=jet
```
![this](img/lowb-t1-aseg.png)

#### Response and fODF estimation

- `dwi2response`
- `dwi2fod`

```bash
dwi2response dwi.mif response.txt -mask mask.mif -sf sf-mask.mif # 36m
dwi2fod dwi.mif response.txt csd.mif -mask mask.mif #
```
This takes about 40 minutes. _TODO images for verification_

If your data was acquired with multiple shells, the above steps should have
warned that the outer shell was automatically selected. We can then estimate
response per shell
```bash
for b in $(mrinfo -shells dwi.mif)
do
  amp2sh dwi.mif dwiSH-$b.mif -shell $b
  maskfilter mask.mif erode - -npass 3 | dwi2tensor dwi.mif tensor-$b.mif -mask -
  tensor2metric tensor-$b.mif -fa - | mrthreshold - sf-$b.mif -abs 0.7
  tensor2metric tensor-$b.mif -vector dirs-$b.mif
  sh2response dwiSH-$b.mif sf-$b.mif dirs-$b.mif response-$b.txt
done
```
_to be tested!_

*TODO* Look at [global tractography](https://github.com/MRtrix3/mrtrix3/wiki/Global-Tractography).
Requires estimating multishell estimations, and use of the `updated_syntax` branch.

#### Anatomical constraints

MRtrix uses the grey matter - white matter interface to improve
track seeding and so-called ACT with several criteria on tissue
type to improve track termination, so generate the necessary volumes:
```bash
act_anat_prepare_fsl t1-in-d.nii.gz act.mif
5tt2gmwmi act.mif gmwmi.mif
```
These steps take 9 minutes.

### Connectivity

The previous section works with volumetric diffusion weighted images. Obtaining
connectomes from such images requires generating tracks and applying volumetric
parcellations to tracks.

#### Tractography

```bash
opt="-seed_gmwmi gmwmi.mif -act act.mif -unidirectional -maxlength 250 -step 0.5"
tckgen $opt -num 10000000 csd.mif brain.tck
```
About 200 tracks per second are selected, so a 10M track tractogram
requires about 13 hours. Again, this can be accelerated by allocated more
CPUs to your VM and changing the MRtrix configuration as described above.

_This is one step which is significantly slower in a VM: natively, on a
Haswell 4 core Xeon, running tckgen with 8 threads, can select ~4.2k
tracks per second, so 10M takes less than an hour._

SIFT filters tracks to improve the fit with the diffusion image,
```bash
tcksift brain.tck csd.mif brain-sift.tck -act act.mif -term_number 2000000
```
This requires X hours.  Note that SIFT requires 1 GB of RAM per 4-8M tracks,
so for the default 10M tracks above, you may need to allocate 3 GB of RAM to
your VM. This results in a 2M track tractogram.

#### Generate connectomes

FreeSurfer's strategy for generating parcellations assigns unique ids to different
regions, but the unique list of ids for all regions is not contiguous, so we generate
a secondary parcellation volume with reassigned, contiguous sequence of ids. More
details and example images are provided on the
[relevant MRtrix example wiki page](https://github.com/MRtrix3/mrtrix3/wiki/labelconfig-worked-example).

The file mapping FreeSurfer names & indices to connectome indices is the
_connectome configuration file_.

For the default parcellation, MRtrix ships with such a file, so we can simply
```bash
conf=/usr/local/mrtrix3/src/connectome/config/fs_default.txt
labelconfig aparc+aseg-in-d.nii.gz $conf aparc+aseg-in-d.mif \
  -lut_freesurfer ${FREESURFER_HOME}/FreeSurferColorLUT.txt
```
For the example subdivided, we need to generate a new connectome configuration.
```bash
_TODO Python script for this_. Then, run `labelconfig` for each custom parcellation,
```
for parc in aparc250
do
  for hemi in lh rh
  do
python annot_to_lut_conn_conf.py lh aparc250 aparc250-lut.txt aparc250-conf.txt
labelconfig $parc+aseg-in-d.nii.gz $conf $parc+aseg-in-d.mif \
  -lut_freesurfer ${FREESURFER_HOME}/FreeSurferColorLUT.txt
done
```

Once the connectome-configured parcellation volumes are prepared, we can
apply them to the tractogram generated above, obtained matrices of
track counts and mean track lengths, for each parcellation volume:
```bash
assignment="-assignment_radial_search 2"
for parc in aparc aparc250
do
  for metric in count meanlength
  do
    tck2connectome brain-sift.tck $parc+aseg-in-d.mif \
      $parc-$metric.csv $assignment -zero_diagonal -metric $metric
  done
done
```

#### Track length distributions

One curiosity of tractography is that regions may have quite wide
track length distributions, which may have an effect on simulation of
conduction velocity. These distributions can be obtained by dumping
lengths & node assignments,

```bash
tckstats brain.tck  -dump lengths.txt
tck2connectome brain-sifted.tck aparc+aseg_d.mif weights.txt -out_assignments assignments.txt
```
and analyzing with a custom script
```bash
python analyze_track_lengths.py
```
