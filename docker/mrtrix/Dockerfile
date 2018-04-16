FROM debian:9

#RUN yum update && yum install wget -y

#RUN apt-get update && apt-get -y install centos-release-scl

RUN apt-get update && apt-get install -y git g++ python python-numpy libeigen3-dev zlib1g-dev libqt4-opengl-dev libgl1-mesa-dev libfftw3-dev libtiff5-dev
#RUN yum install -y devtoolset-7-gcc*

#RUN export PATH=$PATH:/usr/lib64/qt5/bin/

RUN cd /opt && git clone https://github.com/MRtrix3/mrtrix3.git
RUN cd /opt/mrtrix3 && export EIGEN_CFLAGS="-isystem /usr/include/eigen3" && ./configure
RUN cd /opt/mrtrix3 && NUMBER_OF_PROCESSORS=1 ./build
