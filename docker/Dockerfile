FROM debian:9

RUN apt-get update && apt-get install wget -y
RUN wget --no-verbose ftp://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/dev/freesurfer-linux-centos7_x86_64-dev.tar.gz
RUN tar -xzf freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.0.tar.gz -C /opt

COPY license.txt /opt/freesurfer/license.txt

ENV FREESURFER_HOME /opt/freesurfer
RUN source $FREESURFER_HOME/SetUpFreeSurfer.sh

CMD ['mri_info']
