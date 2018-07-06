#$(condor_config_val MASTER) -f -t >> /var/log/condor/MasterLog 2>&1
echo "123456" | sudo -S condor_master
sleep 30

cd pegasus
python run_sequentially.py $1