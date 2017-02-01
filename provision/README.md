# Setup

This requires Vagrant & VirtualBox. Trust me, it's gonna be great.

- get freesurfer & mne packages (see below)
- copy them to `provision/` folder
- run `vagrant up`
- go get a coffee, maybe lunch
- see `provision/jupyter-lab.txt` for URL to open

What do it do though?

- sets up / compiles FSL, FreeSurfer, MNE, OpenMEEG, MRtrix3
- & Python 3.6 w/ all required packages
- starts a Jupyter Notebook & Lab instance to interact with software

The software is all available 

## outline

- freesurfer -> bem/headmat, dmr/trac, decim surf, ct/elec
- subparc & trac -> conn
- conn & dceim surf & headmat -> eeg/meg/bold fwd
- & ct/elec -> seeg fwd

## Get FreeSurfer & MNE

For academic reasons, these softwares require registration to get them. Place
the packages in the `provision/` subdirectory.

If you have access to our cluster, you should be able to copy them from the 
`/soft/dl/` folder, but you still need to get a FreeSurfer license.

### Freesurfer

First, download FreeSurfer v6 for Linux

ftp://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/6.0.0/freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.0.tar.gz

This can take a while as FreeSurfer is a large download.

While it's downloading, you should complete registration at 

https://surfer.nmr.mgh.harvard.edu/fswiki/Registration

(which is free), and they will send you a small license file with some text like
```
foo@bar.com
123456
 *jf0a0fj2932j
 jsrl900r9u309
```
which you should place in a text file at `provision/freesurfer-license.txt`.

#### MNE tools

MNE tools are yet another package which requires manual registration and download from

http://www.nmr.mgh.harvard.edu/martinos/userInfo/data/MNE_register/index.php

After registration, enter your email and password, and click `Download`. On the
resulting page, click the `MNE-2.7.0-3106-Linux-x86_64.tar.gz` link.

Once the download is finished, place this file in the `provision` directory

