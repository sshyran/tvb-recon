# Setup

- setting up a VM with Cent OS 7
- scripts to handle openmeeg, mrtrix compilation
- each step with example results & time estimates

- freesurfer -> bem/headmat, dmr/trac, decim surf, ct/elec
- subparc & trac -> conn
- conn & dceim surf & headmat -> eeg/meg/bold fwd
- & ct/elec -> seeg fwd

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

