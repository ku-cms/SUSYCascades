make clean
make cmssw -j 8
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Summer20UL16APV_106X_SMS.list --verbose --sms
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Summer20UL16_106X_SMS.list --verbose --sms
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Summer20UL17_106X_SMS.list --verbose --sms
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Summer20UL18_106X_SMS.list --verbose --sms
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Summer22_130X_SMS.list --verbose --sms --privateMC
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Summer22_130X_Cascades.list --verbose --cascades --privateMC
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Summer23BPix_130X_Cascades.list --verbose --cascades
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Summer20UL16APV_106X.list --verbose
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Summer20UL16_106X.list --verbose
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Summer20UL17_106X.list --verbose
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Summer20UL18_106X.list --verbose
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Summer22_130X.list --verbose
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Summer22EE_130X.list --verbose
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Summer23_130X.list --verbose
python3 scripts/condor_submit_nano_connect_histograms.py -list samples/NANO/Lists/Summer23BPix_130X.list --verbose
