# -*- mode: ruby -*-
# vi: set ft=ruby :

# Vagrant boxes are usually 8 GB, so we create a separate disk to hold
# work files, cf https://gist.github.com/leifg/4713995. This also the
# file we can share to speed up setup.


# TODO use devpi local cache if doing this alot..

work_disk = 'opt_fs.vdi'

Vagrant.configure("2") do |config|
  config.vm.box = "ubuntu/xenial64"

  config.vm.network "forwarded_port", guest: 8888, host: 8888
  # config.vm.network "private_network", ip: "192.168.33.10"
  # config.vm.synced_folder "../data", "/vagrant_data"

  config.vm.provider "virtualbox" do |vb|
    vb.cpus = 4
    vb.memory = "8192"
    unless File.exist?(work_disk)
      size_in_mb = 50 * 1024
      vb.customize ['createhd', '--filename', work_disk, '--size', size_in_mb]
    end
    vb.customize ['storageattach', :id, '--storagectl', 'SCSI', '--port', 4, '--device', 0, '--type', 'hdd', '--medium', work_disk]
  end

  config.vm.provision "shell", path: "provision/00-system.sh"
  config.vm.provision "shell", privileged: false, path: "provision/10-work-env.sh"

end


