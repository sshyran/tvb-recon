FROM debian:9

#Anaconda
RUN apt-get update && apt-get install -y wget bzip2
RUN wget --quiet --output-document=anaconda.sh "https://repo.continuum.io/archive/Anaconda2-5.1.0-Linux-x86_64.sh" && sh anaconda.sh -b -p /opt/conda
#    && rm -r miniconda.sh


#Mrtrix
RUN apt-get update && apt-get install -y git g++ python python-numpy libeigen3-dev zlib1g-dev libqt4-opengl-dev libgl1-mesa-dev libfftw3-dev libtiff5-dev

RUN cd /opt && git clone https://github.com/MRtrix3/mrtrix3.git
RUN cd /opt/mrtrix3 && export EIGEN_CFLAGS="-isystem /usr/include/eigen3" && ./configure
RUN cd /opt/mrtrix3 && NUMBER_OF_PROCESSORS=1 ./build


#Freesurfer
#RUN apt-get update && apt-get install -y wget

RUN wget --no-verbose ftp://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/dev/freesurfer-linux-centos7_x86_64-dev.tar.gz
RUN tar -xzf freesurfer-linux-centos7_x86_64-dev.tar.gz -C /opt

COPY license.txt /opt/freesurfer/license.txt

ENV FREESURFER_HOME /opt/freesurfer
#RUN sh ${FREESURFER_HOME}/SetUpFreeSurfer.sh


#FSL
RUN wget https://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py

RUN echo "" | python fslinstaller.py
ENV FSLDIR /usr/local/fsl
RUN . ${FSLDIR}/etc/fslconf/fsl.sh
ENV PATH ${FSLDIR}/bin:${PATH}


#TVB-RECON
RUN cd /opt && git clone https://github.com/the-virtual-brain/tvb-recon.git

ENV PATH /opt/conda/bin:${PATH}

RUN conda install -y setuptools numpy scipy matplotlib pytest h5py scikit-learn Cython # trimesh anytree gdist

RUN apt-get install -y python-pip

RUN cd /opt/tvb-recon && python setup.py develop

#Pegasus and HTCondor
RUN apt-get update && apt-get install -y gnupg
RUN wget -O - http://download.pegasus.isi.edu/pegasus/gpg.txt | apt-key add -
RUN echo 'deb http://download.pegasus.isi.edu/pegasus/debian stretch main' >/etc/apt/sources.list.d/pegasus.list
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y condor && apt-get install -y pegasus

RUN echo "DISCARD_SESSION_KEYRING_ON_STARTUP=False\n" >> /etc/condor/condor_config.local
RUN echo "USE_PROCD = False\n" >> /etc/condor/condor_config.local

RUN conda create -n tvb_recon_python3_env python=3.6 anaconda
RUN conda install graphviz
RUN apt-get install -y vim sudo

ENV SUBMIT_USER submitter
ENV GID 1000
ENV UID 1000
ENV PASS 123456

RUN groupadd -g ${GID} ${SUBMIT_USER} && \
    useradd -m -u ${UID} -g ${GID} ${SUBMIT_USER} && \
    usermod -a -G condor ${SUBMIT_USER} && \
    usermod -a -G sudo ${SUBMIT_USER} && \
    echo ${SUBMIT_USER}:${PASS} | chpasswd #&& \
RUN sed -i 's/\(^Defaults.*requiretty.*\)/#\1/' /etc/sudoers
RUN chown -R ${SUBMIT_USER}:${SUBMIT_USER} /home/${SUBMIT_USER}
RUN chmod -R 777 /home/${SUBMIT_USER}

RUN mkdir /opt/pegasus-run && mkdir /opt/pegasus-run/submit && mkdir /opt/pegasus-run/scratch
RUN chown -R ${SUBMIT_USER}:${SUBMIT_USER} /opt/pegasus-run

RUN apt-get install -y tcsh bc

RUN chmod -R 777 /opt/freesurfer/subjects


WORKDIR /opt/tvb-recon
USER ${SUBMIT_USER}

