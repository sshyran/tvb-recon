FROM fs_trial

RUN apt-get update && apt-get install -y git g++ python python-numpy libeigen3-dev zlib1g-dev libqt4-opengl-dev libgl1-mesa-dev libfftw3-dev libtiff5-dev

RUN cd /opt && git clone https://github.com/MRtrix3/mrtrix3.git
RUN cd /opt/mrtrix3 && export EIGEN_CFLAGS="-isystem /usr/include/eigen3" && ./configure
RUN cd /opt/mrtrix3 && NUMBER_OF_PROCESSORS=1 ./build

RUN apt-get update && apt-get install -y wget

RUN wget https://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py

RUN echo "" | python fslinstaller.py 
ENV FSLDIR /usr/local/fsl
RUN . ${FSLDIR}/etc/fslconf/fsl.sh
ENV PATH ${FSLDIR}/bin:${PATH}

#RUN apt-get install -y python-setuptools

RUN cd /opt && git clone https://github.com/the-virtual-brain/tvb-recon.git

RUN wget https://repo.continuum.io/archive/Anaconda2-5.1.0-Linux-x86_64.sh

RUN . Anaconda2-5.1.0-Linux-x86_64.sh -b
#ENV PATH ${HOME}/anaconda2/bin:${PATH}

#RUN conda install -y setuptools numpy scipy matplotlib pytest h5py  # python-scikit-learn python-trimesh python-anytree python-Cython python-gdist

#RUN cd /opt/tvb-recon && python setup.py develop
