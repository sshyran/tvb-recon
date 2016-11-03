#!/bin/bash

git clone https://github.com/MRtrix3/mrtrix3
pushd mrtrix3
MOC=moc-qt4 QMAKE=qmake-qt4 ./configure
./build
popd
cp -r mrtrix3 /usr/local/

cat > /etc/mrtrix.conf <<EOF
NumberOfThreads: 1
EOF

