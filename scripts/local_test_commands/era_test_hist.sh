make clean
make cmssw -j 8
# UL 17
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Additional_Lists/Summer20UL17_106X_local.list --verbose --force --single
# UL 18
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Additional_Lists/Summer20UL18_106X_local.list --verbose --force --single
# Run3 23BPix
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Additional_Lists/Summer23BPix_130X_local.list --verbose --force --single
# Run3 Cascades (Central)
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Additional_Lists/Summer23BPix_130X_Cascades_local.list --verbose --cascades --force --single
# Run3 Cascades (Private)
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Additional_Lists/Summer22_130X_Cascades_local.list --verbose --cascades --privateMC --force --single
# Run3 22SMS
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Additional_Lists/Summer22_130X_SMS_local.list --verbose --privateMC --sms --force --single
