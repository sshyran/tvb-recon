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
RUN apt-get install -y tcsh bc libgomp1 perl-modules

RUN mkdir /opt/freesurfer-stable && wget -N -qO- ftp://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/6.0.0/freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.0.tar.gz | tar -xzv -C /opt/freesurfer-stable
COPY license.txt /opt/freesurfer-stable/license.txt



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

COPY license.txt /opt/freesurfer-stable/freesurfer/license.txt
RUN chmod -R 777 /opt/freesurfer-stable/freesurfer/subjects

# Have a generic data zone for mounted volume
#RUN mkdir /home/submitter/Cleveland_datasets && mkdir /home/submitter/Cleveland_datasets/subjects
#RUN chmod -R 777 /home/submitter/Cleveland_datasets/subjects

RUN sudo ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh
RUN chmod -R 777 /opt/freesurfer-stable/freesurfer/bin

#ENV OS Linux
#ENV FS_OVERRIDE 0
#ENV SUBJECTS_DIR /home/submitter/Cleveland_datasets/subjects
#ENV FSF_OUTPUT_FORMAT nii.gz
#ENV MNI_DIR /opt/freesurfer/mni
#ENV LOCAL_DIR /opt/freesurfer/local
#ENV FREESURFER_HOME /opt/freesurfer
#ENV FSFAST_HOME /opt/freesurfer/fsfast
#ENV MINC_BIN_DIR /opt/freesurfer/mni/bin
#ENV MINC_LIB_DIR /opt/freesurfer/mni/lib
#ENV MNI_DATAPATH /opt/freesurfer/mni/data
#ENV FMRI_ANALYSIS_DIR /opt/freesurfer/fsfast
#ENV PERL5LIB /opt/freesurfer/mni/lib/perl5/5.8.5
#ENV MNI_PERL5LIB /opt/freesurfer/mni/lib/perl5/5.8.5
#ENV PATH /opt/freesurfer/bin:/opt/freesurfer/fsfast/bin:/opt/freesurfer/tktools:/opt/freesurfer/mni/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
#RUN /bin/bash -c ' echo -e "source $FREESURFER_HOME/FreeSurferEnv.sh &>/dev/null" >> /root/.bashrc '

#COPY start-condor.sh /usr/sbin
#RUN chmod -R 777 /usr/sbin/start-condor.sh

WORKDIR /opt/tvb-recon

RUN apt-get install -y python-pip
RUN [ "/bin/bash", "-c", "source /opt/conda/bin/activate tvb_recon_python3_env && python setup.py develop && pip install nibabel" ]

RUN apt-get install -y dc

USER ${SUBMIT_USER}


#CMD ["/usr/sbin/start-condor.sh"]

#RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
#RUN conda activate tvb_recon_python3_env && python setup.py develop && conda deactivate
