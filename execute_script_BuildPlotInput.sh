#!/bin/bash
tar -xzf ./config_BuildPlotInput.tgz

export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh

source ./config_BuildPlotInput/cmssw_setup_connect.sh

cmssw_setup sandbox-CMSSW_10_6_5-6403d6f.tar.bz2

source ./config_BuildPlotInput/setup_RestFrames_connect.sh 

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-pafccj3/lib/

./config_BuildPlotInput/BuildPlotInput.x "$@"
