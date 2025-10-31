# runs skim to check that single test condor job works

make clean
make cmssw -j 8

### preUL
###python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Fall17_102X_local.list --verbose --force --single --sys

# preUL SMS
#python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Fall17_102X_SMS_T1_local.list --verbose --sms --force --single --sys

# UL 17
python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Summer20UL17_106X_local.list --verbose --force --single --sys

# UL 18
python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Summer20UL18_106X_local.list --verbose --force --single --sys

# Run3 23BPix
python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Summer23BPix_130X_local.list --verbose --force --single --sys

# Run3 Cascades (Central)
python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Summer23BPix_130X_Cascades_local.list --verbose --cascades --force --single --sys

# Run3 Cascades (Private)
python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Summer22_130X_Cascades_local.list --verbose --cascades --privateMC --force --single --sys

# Run2 Data
python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Summer20UL17_106X_Data_local.list --verbose --data --force --single 

# Run3 Data
python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Summer23BPix_130X_Data_local.list --verbose --data --force --single 

echo "Finished running era test!"
