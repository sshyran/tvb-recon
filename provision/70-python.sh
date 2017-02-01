#!/usr/bin/env bash

# bootstrap recent Python, assuming sane dev env & ssl headers
set -eu
set -o pipefail

export PREFIX=${PREFIX:-"/work/env"}

if [[ ! -z $(which jupyter) ]]
then
    echo "[70-python.sh] jupyter found, not setting up Python env."
    exit 0
else
    echo "[70-python.sh] building Python env."
fi

j=6
prefix=/work/env

zver=1.2.11
zlib_url=http://zlib.net/zlib-$zver.tar.gz

bzver=1.0.6
bzip_url=http://www.bzip.org/$bzver/bzip2-$bzver.tar.gz

sql_url=https://sqlite.org/2017/sqlite-autoconf-3160200.tar.gz

pyver=3.6.0
Py_url=https://www.python.org/ftp/python/$pyver/Python-$pyver.tgz

for pkg in zlib bzip sql Py
do
    curl -O $(eval echo \$${pkg}_url)
    tar xzf ${pkg}*
    pushd ${pkg}*
    if [[ $pkg == "bzip" ]]
    then
        make -j$j install PREFIX=$prefix
    else
        ./configure --prefix=$prefix
        make -j$j
        make install
    fi
    popd
done


py_pkgs="numpy scipy matplotlib cython scikit-learn pandas h5py nibabel"
py_pkgs="$py_pkgs nibabel anytree trimesh flake8 mypy mne jupyterlab"
for pkg in $py_pkgs
do
    echo "pip installing $pkg"
    $prefix/bin/pip3 install $pkg
done

$prefix/bin/jupyter serverextension enable --py jupyterlab --sys-prefix

# consider focusing on jupyter/notebook for visualizations are required

# finally set up recon tools package for work
echo TODO pushd /vagrant
echo TODO $prefix/bin/python setup.py develop
