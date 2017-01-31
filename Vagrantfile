# -*- mode: ruby -*-
# vi: set ft=ruby :

# Vagrant boxes are usually 8 GB, so we create a separate disk to hold
# work files, cf https://gist.github.com/leifg/4713995. This also the
# file we can share to speed up setup.
work_disk = 'opt_fs.vdi'

Vagrant.configure("2") do |config|
  config.vm.box = "ubuntu/xenial64"

  config.vm.network "forwarded_port", guest: 8888, host: 8888
  # config.vm.network "private_network", ip: "192.168.33.10"
  # config.vm.synced_folder "../data", "/vagrant_data"


  # TODO ansible this / make reusable for cluster
  # TODO CI/CD the setup, make provisioning faster
  config.vm.provider "virtualbox" do |vb|
    vb.cpus = 4
    vb.memory = "8192"
    unless File.exist?(work_disk)
      size_in_mb = 50 * 1024
      vb.customize ['createhd', '--filename', work_disk, '--size', size_in_mb]
    end
    vb.customize ['storageattach', :id, '--storagectl', 'SCSI', '--port', 4, '--device', 0, '--type', 'hdd', '--medium', work_disk]
  end

  config.vm.provision "shell", inline: <<-SHELL
    apt-get update && apt-get upgrade -y
    apt-get install -y tcsh libeigen3-dev liblapack-dev libblas-dev libssl-dev

    # setup partition
    parted /dev/sdc mklabel msdos
    parted /dev/sdc mkpart primary 512 100%
    mkfs.xfs /dev/sdc1
    mkdir /work
    echo "/dev/sdc1 /work xfs noatime,nobarrier 0 0" >> /etc/fstab
    mount /work
    mkdir /work/{data,env}
    chown -R ubuntu:ubuntu /work

    # setup /work/env/lib as system wide library location
    echo /work/env/lib > /etc/ld.so.conf.d/work.conf
    ldconfig

    cp /vagrant/jupyter-lab.service /etc/systemd/system/

    # system-wide bash env
    # TODO debug this
    echo 'export PREFIX=/work/env' >> /etc/bash.bashrc
    echo 'export FREESURFER_HOME=$PREFIX/freesurfer' >> /etc/bash.bashrc
    echo 'export SUBJECTS_DIR=/work/data' >> /etc/bash.bashrc
    echo 'export PATH=$PREFIX/bin:$PATH' >> /etc/bash.bashrc
    echo 'source $FREESURFER_HOME/FreeSurferEnv.sh' >> /etc/bash.bashrc
  SHELL

  config.vm.provision "shell", privileged: false, inline: <<-SHELL
    exit 0
    pushd /work/env

    bash /vagrant/docs/setup-vm/70-python.sh
    bash /vagrant/docs/setup-vm/10-cmake.sh
    bash /vagrant/docs/setup-vm/35-openmeeg.sh

    tar -C /work/env -xzf /vagrant/archives/freesurfer*.tar.gz

    git clone https://github.com/mrtrix3/mrtrix3 /work/env/mrtrix3
    pushd /work/env/mrtrix3
    ./configure -nogui
    ./build
    popd

    sudo systemctl start jupyter-lab
    sleep 5
    sudo journalctl -u jupyter-lab

    # TODO MNE
    # TODO check all OK
  SHELL

end



