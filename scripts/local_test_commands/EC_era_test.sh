# runs EC locally to check that EC things work locally

make clean
make cmssw -j 8

# Run3 Cascades (Central)
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Additional_Lists/Summer23BPix_130X_Cascades_local.list --cascades --dryrun 
./MakeEventCount_NANO.x -ilist=./Summer23BPix_130X_Cascades_local_EventCount/config/list/SlepSnuCascade_MN1-220_MN2-260_MC1-240_TuneCP5_13p6TeV_madgraphMLM-pythia8_Summer23BPix_130X/SlepSnuCascade_MN1-220_MN2-260_MC1-240_TuneCP5_13p6TeV_madgraphMLM-pythia8_Summer23BPix_130X_0.list -ofile=EC_SlepSnuCascade_MN1-220_MN2-260_MC1-240_TuneCP5_13p6TeV_madgraphMLM-pythia8_Summer23BPix_130X_0.root -tree=Events --cascades -dataset=SlepSnuCascade_MN1-220_MN2-260_MC1-240_TuneCP5_13p6TeV_madgraphMLM-pythia8 -filetag=Summer23BPix_130X

# Run3 Cascades (Private)
python3 scripts/condor_submit_nano_cmslpc_eventCount.py --connect -list samples/NANO/Lists/Additional_Lists/Summer22_130X_Cascades_local.list --cascades --privateMC --dryrun
./MakeEventCount_NANO.x -ilist=./Summer22_130X_Cascades_local_EventCount/config/list/SlepSnuCascade_220-209_200-190-180_2022_NANO_JustinPrivateMC_Summer22_130X_Cascades_Summer22_130X/SlepSnuCascade_220-209_200-190-180_2022_NANO_JustinPrivateMC_Summer22_130X_Cascades_Summer22_130X_0.list -ofile=EC_SlepSnuCascade_220-209_200-190-180_2022_NANO_JustinPrivateMC_Summer22_130X_Cascades_Summer22_130X_0.root -tree=Events --cascades --privateMC -dataset=SlepSnuCascade_220-209_200-190-180_2022_NANO_JustinPrivateMC_Summer22_130X_Cascades -filetag=Summer22_130X 

echo "Finished running EC era test!"
