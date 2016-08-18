#!/bin/bash

for script in 00-system 10-cmake
do
	bash ${TVBVIRT}/setup/${script}.sh
done

echo "rebooting in 5 seconds.."
sleep 5
reboot
