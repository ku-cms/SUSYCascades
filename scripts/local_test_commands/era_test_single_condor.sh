# runs skim locally to check that things work locally

make clean
make cmssw -j 8

### preUL
###python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Fall17_102X_local.list --verbose --force --single

# preUL SMS
python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Fall17_102X_SMS_T1_local.list --verbose --sms --force --single

# UL 17
python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Summer20UL17_106X_local.list --verbose --force --single

# UL 18
python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Summer20UL18_106X_local.list --verbose --force --single

# Run3 23BPix
python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Summer23BPix_130X_local.list --verbose --force --single

# Run3 Cascades (Central)
python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Summer23BPix_130X_Cascades_local.list --verbose --cascades --force --single

# Run3 Cascades (Private)
python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Summer22_130X_Cascades_local.list --verbose --cascades --privateMC --force --single

#### Use this set of options to process central cascades MC like old way (not recommended)
###python3 scripts/condor_submit_nano_connect_ntuples.py -list samples/NANO/Lists/Additional_Lists/Summer23BPix_130X_SMS_local.list --verbose --sms --force --single

echo "Finished running era test!"
