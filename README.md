# SUSYCascades
**SUSYCascades** code package implementing the KU Cascades Analysis

---------------------
CMSSW-dependent build 
---------------------

In order for the SUSYCascades package to build the `BuildFit.x` executable,
the user must make sure that **CMSSW** is available and include the **SUSYCascades** package
in the correct location in the **CMSSW** directory structure.
<!--- You will also need the **CombineHarvester** and **HiggsAnalysis** **CMSSW** packages.
These packages must be included in the **CMSSW** directory structure as:
```
- CMSSW_Z_Y_X
  - src
    - CombineHarvester
    - HiggsAnalysis
    - SUSYCascades
 ```
--->
You can setup a **CMSSW** area and checkout the required packages by performing the terminal commands below. 
These instructions are for bash users.
Bash is recommended, as framework scripts for other shells are not created/supported.
You can check what shell you are using by running `echo $0` or `echo $SHELL`.

### Setup environmental variables for CMSSW
Put these commands in either `~/.bash_profile` or `~/.bashrc` (but not both)
so that they are run at the beginning of every new terminal session:
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/cms.cern.ch/crab3/crab.sh
```    
Once these commands are added to `~/.bash_profile` or `~/.bashrc` for the first time,
you can then source whichever file was modified to apply these changes.
```
source ~/.bash_profile
source ~/.bashrc
```

### Setup CMSSW area
Checkout CMSSW_13_3_1.
```
cmsrel CMSSW_13_3_1
cd CMSSW_13_3_1/src
cmsenv
```
    
### Checkout required packages
Download this repos in the `CMSSW_13_3_1/src` directory.
<!--- Use the KU branch of HiggsAnalysis, the connect branch of CombineHarvester, and the master branch of SUSYCascades.
```
 git clone -b KU https://github.com/zflowers/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
 git clone -b connect https://github.com/zflowers/CombineHarvester.git CombineHarvester
 git clone -b master https://github.com/ku-cms/SUSYCascades.git SUSYCascades
```
--->
```
 git clone -b master https://github.com/ku-cms/SUSYCascades.git SUSYCascades
```

### FNAL CMS LPC
If running on FNAL cmslpc, then the command to copy the tar ball is:
```
cp /uscms/home/z374f439/nobackup/RestFrames_vNewest.tar ./
```

### CMS Connect
If running on CMS connect, then the command to copy the tar ball is:
```
cp /ospool/cms-user/zflowers/public/RestFrames_vNewest.tar ./
```

### LXPLUS
If running on lxplus, then the command to copy the tar ball is:
```
cp /afs/cern.ch/user/z/zflowers/storage/public/RestFrames_vNewest.tar ./
```

### Setup RestFrames
Setup the RestFrames framework in the `CMSSW_13_3_1/src` directory.
You'll copy the tar ball from the area below and unpack it
```
tar -xvf RestFrames_vNewest.tar
cd RestFrames
./configure --prefix=$CMSSW_BASE/src/restframe_build
make -j8
make install
cd ..
```

### Build/compile
Build the CMSSW packages from the `CMSSW_13_3_1/src` directory.
```
scram b -j8
```
Then build **SUSYCascades** (with the `BuildFit.x` executable).
Note that you will need to run the RestFrames setup script
(in this case the one specific to CMS connect) before compiling.
Also note that we compile using CMSSW, which is required.
Compile and run issues have been observed when compiling without CMSSW.
```
cd SUSYCascades
source scripts/setup_RestFrames_connect.sh
make clean && make cmssw -j8
```
### To submit ntuples:
Choose appropriate list from samples/NANO/Lists/
If there are datasets that you do not want to process, you can comment out with '#'
To submit use scripts/condor_submit_nano_connect_ntuples.py
```
python3 scripts/condor_submit_nano_connect_ntuples.py -split 20 -list samples/NANO/Lists/NAME_OF_YOUR_LIST.list
```
If running on data or signal, add --sms or --data
-split N defines how many jobs to split each input file into.
--count can be useful before launching jobs to check how many jobs you will be submitting.
NOTE: max number of jobs per batch (dataset) is ~15k, max number of jobs per user is 100k.
The submission script will protect against going over the second limit (user total).
--csv will add a csv into the submission working directory will job info (not needed but useful).

### After NTUPLE jobs finish
Important to check for job failures after jobs finish.
Note that you do not have to wait until EVERY job finishes before starting to run checks.
At the end of the submit script there is a print statement that tells you how to run the checker (also posted below)
```
nohup python3 python/CheckFiles.py -d f{filetag}/ -o {OUT_DIR} -e > CheckFiles_{filetag}.debug 2>&1 &
```
Recommended to run in background with nohup and direct output to a log file.
For ~100 datasets can take up to 8 hours the first time if all checks are run.
Use python3 python/CheckFiles.py to see breakdown of args.
-e is recommended because it is a) slow and b) can return false negatives.

### After NTUPLE jobs pass checks
After ntuple jobs are finished, should use hadd to reduce number of files.
Python script in scripts/DO_hadd.py can be used to efficiently combine output files.
Has guardrails in place to handle large number of output files (another reason to keep total jobs/batch under 15k) and larger sets of data (>100 GB)
Example submission:
```
nohup python3 scripts/DO_hadd.py -idir /ospool/cms-user/$USER/NTUPLES/Processing/{filetag}/ -odir /local-scratch/$USER/NTUPLES/HADD/{filetag}/ > HADD_logs/HADD_{filetag}.debug 2>&1 &
```
Again, note the use of nohup, &, and >.

### Making EventCount files
EventCount files (ECs) store the total number of events and sum of weights for each event in every sample.
They are needed in order to properly calculate weights in reduced ntuples.
There are a few python scripts intended for batch updating all samples from the main lists (needed for both bkg and signal).
Run this to batch submit EC jobs
```
nohup ./scripts/submit_event_count.sh > submit_EC.debug 2>&1 &
```
After EC jobs finish. Need to check for failed jobs.
```
nohup python3 python/DASCheck_EventCount.py -d > /dev/null 2>&1 &
```
Finally run HADD script to update final EC files used by ntuple jobs.
```
nohup python3 scripts/HADD_EventCount.py > HADD_EventCount.debug 2>&1 &
```

#### User quotas
Connect:

/home/ = 100GB

/ospool/ = 1TB

/local-scratch/ = 1TB

### Updating Lists
If needed to update lists or update cross-sections use this repo: https://github.com/ku-cms/ListMaker
Only use if dataset is not already in samples/NANO/, otherwise can modify an existing list

<!---
 ### Running combineTool.py on cms connect syntax
 To run text2workspace (T2W):
 Go to the directory with the datacard for a given mass point
 ```
 combineTool.py -M T2W -i datacard.txt -m MASS_Value --job-mode connect -o workspace.root --input-file ../../../FitInput_KUEWKino_2017.root --sub-opts='+ProjectName=cms.org.ku" \n request_memory = 8 GB \n'
 ```
 
 Note that mass value will typically be the name of the directory that the datacard is in (ex: 5000325) and the input file is the path to the output root file from BuildFit.x
 
 To run AsymptoticLimits:
 Go to the directory that was the output of BuildFit.x 
 ```
 combineTool.py -M AsymptoticLimits -d */*/*/datacard.txt --job-mode connect --input-file FitInput_KUEWKino_2017.root --sub-opts='+ProjectName="cms.org.ku" \n request_memory = 8 GB \n'
 ```
 
 To run impacts:
 Go to the directory that the workspace (and datacard) live in
 ```
 combineTool.py -M Impacts -d workspace.root -m MASS_Value --doInitialFit --robustFit 1 --job-mode connect --sub-opts='+ProjectName="cms.org.ku" \n request_memory = 8 GB \n'
 combineTool.py -M Impacts -d workspace.root -m 5000325 --input-file ../ --job-mode connect --sub-opts='+ProjectName="cms.org.ku" \n request_memory = 8 GB \n' 
 ```
 
 Note that the second step only should be ran after the jobs in the first step finish
 The input file is just the entire directory that the workspace lives in
--->

