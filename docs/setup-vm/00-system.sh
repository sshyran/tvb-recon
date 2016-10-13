#!/bin/bash

yum update -y
yum install -y gcc gcc-c++ cmake git blas-devel atlas-devel zlib-devel bzip2 tcsh numpy python-matplotlib scipy python-setuptools bc Cython mesa-libGLU libXScrnSaver qt qt-devel
easy_install nibabel
easy_install gdist

echo /usr/local/lib >> /etc/ld.so.conf
echo /usr/local/lib64 >> /etc/ld.so.conf
ldconfig

yum group install -y "X Window System"
yum install -y gnome-classic-session gnome-terminal nautilus-open-terminal control-center liberation-mono-fonts firefox libpng12 libmng
ln -s /usr/lib64/libpangoxft-1.0.so.0 /usr/lib64/libpangox-1.0.so.0
unlink /etc/systemd/system/default.target
ln -sf /lib/systemd/system/graphical.target /etc/systemd/system/default.target
