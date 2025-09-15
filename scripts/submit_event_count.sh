# Use this script to submit ALL EC jobs
make clean
make cmssw -j 8
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer16_102X.list
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer16_102X_SMS.list --sms
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Fall17_102X.list
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Fall17_102X_SMS.list --sms
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Autumn18_102X.list
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Autumn18_102X_SMS.list --sms
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer20UL16_106X.list
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer20UL16_106X_SMS.list --sms
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer20UL16APV_106X.list
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer20UL16APV_106X_SMS.list --sms
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer20UL17_106X.list
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer20UL17_106X_SMS.list --sms
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer20UL18_106X.list
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer20UL18_106X_SMS.list --sms
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer22_130X.list
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer22EE_130X.list
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer23_130X.list
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer23BPix_130X.list
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer23BPix_130X_SMS.list --sms --cascades
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer23BPix_130X_Cascades.list --cascades
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer22_130X_Cascades.list --cascades --privateMC
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Summer22_130X_SMS.list --sms --privateMC
# Run checker after submitting all jobs
nohup python3 python/DASCheck_EventCount.py -d > /dev/null 2>&1 &
