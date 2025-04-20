source /cvmfs/cms.cern.ch/cmsset_default.sh
tar -xzf ./config.tgz
source ./config/cmssw_setup_connect_el9.sh
cmssw_setup sandbox-CMSSW_13_3_1-el9.tar.bz2  
source ./config/setup_RestFrames_connect.sh 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/lhapdf/6.4.0-52852f9a177b8e8b5b72e2ae6b1327b6//lib/
export X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates/
echo "ENV setup complete! Running exe..."
./config/MakeReducedNtuple_NANO.x "$@"
