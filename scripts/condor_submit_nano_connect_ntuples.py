#! /usr/bin/env python

import os, sys, time
from colorama import Fore, Back, Style
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))
from CondorJobCountMonitor import CondorJobCountMonitor

# Example submission: 
#  python3 scripts/condor_submit_nano_connect_ntuples.py -split 10 -list samples/NANO/Lists/Fall17_102X.list --sys --slim --csv 
#
#  python3 scripts/condor_submit_nano_connect_ntuples.py -split 20 -list samples/NANO/Lists/Fall17_102X.list --sys --slim --verbose --csv && python3 scripts/condor_submit_nano_connect_ntuples.py -split 30 -list samples/NANO/Lists/Fall17_102X_SMS_FastSim_Quick.list --sys --slim --sms --fastsim --verbose --csv && python3 scripts/condor_submit_nano_connect_ntuples.py -split 10 -list samples/NANO/Lists/Fall17_102X_Data_MET.list --slim --data --verbose --csv 

# ----------------------------------------------------------- #
# Parameters
# ----------------------------------------------------------- #
# current working directory
pwd         = os.environ['PWD']
RUN_DIR     = pwd
jobEXE      = "execute_script.sh"
EXE         = "MakeReducedNtuple_NANO.x"
RESTFRAMES  = './scripts/setup_RestFrames_connect.sh'
#CMSSW_SETUP = './scripts/cmssw_setup_connect.sh'
CMSSW_SETUP = './scripts/cmssw_setup_connect_el9.sh'
TREE        = "Events"
USER        = os.environ['USER']
OUT_BASE    = "/ospool/cms-user/"+USER+"/NTUPLES/Processing"
#OUT_BASE    = "/uscms/home/"+USER+"/nobackup/NTUPLES/Processing"
LIST        = "default.list"
QUEUE       = ""
SPLIT       = 1
THRESHOLD   = 80000
# ----------------------------------------------------------- #

def new_listfile(rootlist, listfile):
    mylist = open(listfile,'w')
    for f in rootlist:
        mylist.write(f+" \n")
    mylist.close()

def create_filelist(rootlist, dataset, filetag):
    listlist = []
    listcount = 0
    
    sublist = []
    for f in rootlist:
        sublist.append(f)
        if len(sublist) >= 1:
            listfile = "%s/%s_%s_%d.list" % (listdir_sam, dataset, filetag, listcount)
            new_listfile(sublist, listfile)
            listlist.append(listfile)
            sublist = []
            listcount += 1

    if len(sublist) > 0:
        listfile = "%s/%s_%s_%d.list" % (listdir_sam, dataset, filetag, listcount)
        new_listfile(sublist, listfile)
        listlist.append(listfile)

    return listlist

def write_sh_single(srcfile,ifile,ofile,logfile,outfile,errfile,dataset,filetag,n,NAME):
    srcfile = srcfile.replace('.submit','_single.submit')
    ofile = ofile.replace('_$(ItemIndex)_$(Step)','_0_0')
    outfile = outfile.replace('_$(ItemIndex)_$(Step)','_0_0')
    errfile = errfile.replace('_$(ItemIndex)_$(Step)','_0_0')
    logfile = logfile.replace('_$(ItemIndex)_$(Step)','_0_0')
    ifile = ifile.replace('_list','_0')
    if DO_SMS == 1:
        ifile = ifile.replace(pwd+'/'+filetag+'_SMS','./config')
    elif DO_DATA == 1:
        ifile = ifile.replace(pwd+'/'+filetag+'_Data','./config')
    else:
        ifile = ifile.replace(pwd+'/'+filetag,'./config')

    fsrc = open(srcfile,'w')
    fsrc.write('# Note: For only submitting 1 job! \n')
    fsrc.write('# output (x_y), list and split args need to be updated \n')
    fsrc.write('# Split should be the second number +1 \n')
    fsrc.write('#       Example: file_1_15 would have a split value of 16,'+str(n)+'\n')
    fsrc.write('universe = vanilla \n')
    fsrc.write('executable = '+jobEXE+" \n")
    fsrc.write('use_x509userproxy = true \n')
    fsrc.write('Arguments = ');
    fsrc.write('-ilist='+ifile+' ')
    fsrc.write('-ofile='+ofile.split('/')[-1]+" ")
    fsrc.write('-tree='+TREE+" ")
    if DO_SMS == 1:
        fsrc.write('--sms ')
    if DO_DATA == 1:
        fsrc.write('--data ')
    if SYS == 1 and DO_DATA != 1:
        fsrc.write('--sys ')
    if FASTSIM == 1 and DO_DATA != 1: # note that technically FS should only be needed for SMS but not requiring it here
        fsrc.write('--fastsim ')
    if SLIM == 1:
        fsrc.write('--slim ')
    fsrc.write('-dataset='+dataset+" ")
    fsrc.write('-filetag='+filetag+" ")
    fsrc.write('-eventcount='+EVTCNT+" ")
    fsrc.write('-filtereff='+FILTEREFF+" ")
    fsrc.write('-json='+JSON+" ")
    fsrc.write('-pu='+PUFOLD+" ")
    fsrc.write('-btag='+BTAGFOLD+" ")
    fsrc.write('-lep='+LEPFOLD+" ")
    fsrc.write('-jme='+JMEFOLD+" ")
    fsrc.write('-metfile='+METFILE+" ")
    fsrc.write('-prefirefile='+PREFIREFILE+" ")
    fsrc.write('-xsjsonfile='+XSJSONFILE+" ")
    fsrc.write('-split=1,'+str(n)+'\n')

    outlog = outfile+".out"
    errlog = errfile+".err"
    loglog = logfile+".log"
    fsrc.write('output = '+outlog+" \n")
    fsrc.write('error = '+errlog+" \n")
    fsrc.write('log = '+loglog+" \n")
    #fsrc.write('Requirements = (Machine != "red-node000.unl.edu" && Machine != "ncm*.hpc.itc.rwth-aachen.de" && Machine != "*mh-epyc7662-8.t2.ucsd.edu" && Machine != "*sdsc-88.t2.ucsd.edu" && Machine != "*beowulf.cluster" && Machine != "*126.hep.olemiss.edu")\n')
    fsrc.write('request_memory = 2 GB \n')
    fsrc.write('+RequiresCVMFS = True \n')
    fsrc.write('priority = 10 \n')
    if USE_URL:
        # Warning: The stash.osgconnect.net endpoint has been decommissioned.
        # CMS connect is working on implementing an OSDF endpoint solution.
        transfer_input = 'transfer_input_files = https://stash.osgconnect.net/cms-user/'+USER+"/"+NAME+"/"+'config.tgz,https://stash.osgconnect.net/cms-user/zflowers/public/sandbox-CMSSW_10_6_5-6403d6f.tar.bz2\n'
    else:
        #transfer_input = 'transfer_input_files = '+TARGET+'config.tgz,/ospool/cms-user/zflowers/public/sandbox-CMSSW_10_6_5-6403d6f.tar.bz2\n'
        transfer_input = 'transfer_input_files = '+TARGET+'config.tgz,/ospool/cms-user/zflowers/public/sandbox-CMSSW_13_3_1-el9.tar.bz2\n'
    fsrc.write(transfer_input)

    fsrc.write('should_transfer_files = YES\n')
    fsrc.write('when_to_transfer_output = ON_EXIT\n')

    transfer_out_files = 'transfer_output_files = '+ofile.split('/')[-1]+'\n'
    fsrc.write(transfer_out_files)

    transfer_out_remap = 'transfer_output_remaps = "'+ofile.split('/')[-1]+'='+ofile
    transfer_out_remap += '"\n'
    fsrc.write(transfer_out_remap)
    # Trying new logic
    fsrc.write('RequestCpus = 1\n')
    fsrc.write('periodic_hold_subcode = 42\n')
    fsrc.write('periodic_hold = (CpusUsage > RequestCpus) && (JobStatus == 2) && (CurrentTime - EnteredCurrentStatus > 60)\n')
    fsrc.write('+JobTransforms = "if HoldReasonSubCode == 42 set RequestCpus = MIN(RequestCpus + 1, 32)"\n')
    fsrc.write('periodic_release = (HoldReasonCode == 12 && HoldReasonSubCode == 256 || HoldReasonCode == 13 && HoldReasonSubCode == 2 || HoldReasonCode == 12 && HoldReasonSubCode == 2 || HoldReasonCode == 26 && HoldReasonSubCode == 120 || HoldReasonCode == 3 && HoldReasonSubCode == 0 || HoldReasonSubCode == 42)\n')
    
    #fsrc.write('RequestCpus=IfThenElse((IsUndefined(CpusUsage) || CpusUsage <= 1),1,MIN(CpusUsage+2, 32))\n')
    #fsrc.write('periodic_hold = (CpusUsage >= RequestCpus) && (JobStatus == 3)\n')
    #fsrc.write('periodic_hold_subcode = 42\n')
    #fsrc.write('periodic_release = (HoldReasonCode == 12 && HoldReasonSubCode == 256 || HoldReasonCode == 13 && HoldReasonSubCode == 2 || HoldReasonCode == 12 && HoldReasonSubCode == 2 || HoldReasonCode == 26 && HoldReasonSubCode == 120 || HoldReasonCode == 3 && HoldReasonSubCode == 0 || HoldReasonSubCode == 42)\n')
    fsrc.write('priority = 10 \n')
    fsrc.write('+REQUIRED_OS="rhel9"\n')
    fsrc.write('job_lease_duration = 3600\n')
    fsrc.write('+ProjectName="cms.org.ku"\n')
    fsrc.write('Requirements = HAS_SINGULARITY == True\n')
    fsrc.write('MY.SingularityImage = "/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel9"\n')
    fsrc.write('queue')
    fsrc.close()

def write_sh(srcfile,ifile,ofile,logfile,outfile,errfile,dataset,filetag,n,NAME):
    fsrc = open(srcfile,'w')
    fsrc.write('universe = vanilla \n')
    fsrc.write('executable = '+jobEXE+" \n")
    fsrc.write('use_x509userproxy = true \n')
    fsrc.write('Arguments = ');
    fsrc.write('-ilist=$(Item) ')
    fsrc.write('-ofile='+ofile.split('/')[-1]+" ")
    fsrc.write('-tree='+TREE+" ")
    if DO_SMS == 1:
        fsrc.write('--sms ')
    if DO_DATA == 1:
        fsrc.write('--data ')
    if SYS == 1 and DO_DATA != 1:
        fsrc.write('--sys ')
    if FASTSIM == 1 and DO_DATA != 1: # note that technically FS should only be needed for SMS but not requiring it here
        fsrc.write('--fastsim ')
    if SLIM == 1:
        fsrc.write('--slim ')
    fsrc.write('-dataset='+dataset+" ")
    fsrc.write('-filetag='+filetag+" ")
    fsrc.write('-eventcount='+EVTCNT+" ")
    fsrc.write('-filtereff='+FILTEREFF+" ")
    fsrc.write('-json='+JSON+" ")
    fsrc.write('-pu='+PUFOLD+" ")
    fsrc.write('-btag='+BTAGFOLD+" ")
    fsrc.write('-lep='+LEPFOLD+" ")
    fsrc.write('-jme='+JMEFOLD+" ")
    fsrc.write('-metfile='+METFILE+" ")
    fsrc.write('-prefirefile='+PREFIREFILE+" ")
    fsrc.write('-xsjsonfile='+XSJSONFILE+" ")
    splitstring = '-split=%s,%d\n' % ('$$([$(Step)+1])', n)
    fsrc.write(splitstring)

    outlog = outfile+".out"
    errlog = errfile+".err"
    loglog = logfile+".log"
    fsrc.write('output = '+outlog+" \n")
    fsrc.write('error = '+errlog+" \n")
    fsrc.write('log = '+loglog+" \n")
    fsrc.write('request_memory = 2 GB \n')
    if USE_URL:
        # Warning: The stash.osgconnect.net endpoint has been decommissioned.
        # CMS connect is working on implementing an OSDF endpoint solution.
        transfer_input = 'transfer_input_files = https://stash.osgconnect.net/cms-user/'+USER+"/"+NAME+"/"+'config.tgz,https://stash.osgconnect.net/cms-user/zflowers/public/sandbox-CMSSW_10_6_5-6403d6f.tar.bz2\n'
    else:
        #transfer_input = 'transfer_input_files = '+TARGET+'config.tgz,/ospool/cms-user/zflowers/public/sandbox-CMSSW_10_6_5-6403d6f.tar.bz2\n'
        transfer_input = 'transfer_input_files = '+TARGET+'config.tgz,/ospool/cms-user/zflowers/public/sandbox-CMSSW_13_3_1-el9.tar.bz2\n'
    fsrc.write(transfer_input)

    fsrc.write('should_transfer_files = YES\n')
    fsrc.write('when_to_transfer_output = ON_EXIT\n')

    transfer_out_files = 'transfer_output_files = '+ofile.split('/')[-1]+'\n'
    fsrc.write(transfer_out_files)

    transfer_out_remap = 'transfer_output_remaps = "'+ofile.split('/')[-1]+'='+ofile
    transfer_out_remap += '"\n'
    fsrc.write(transfer_out_remap)
    
    fsrc.write('+ProjectName="cms.org.ku"\n')
    fsrc.write('+REQUIRED_OS="rhel9"\n')
    fsrc.write('job_lease_duration = 3600\n')
    # Trying new logic
    fsrc.write('RequestCpus = 1\n')
    fsrc.write('periodic_hold_subcode = 42\n')
    fsrc.write('periodic_hold = (CpusUsage > RequestCpus) && (JobStatus == 2) && (CurrentTime - EnteredCurrentStatus > 60)\n')
    fsrc.write('+JobTransforms = "if HoldReasonSubCode == 42 set RequestCpus = MIN(RequestCpus + 1, 32)"\n')
    fsrc.write('periodic_release = (HoldReasonCode == 12 && HoldReasonSubCode == 256 || HoldReasonCode == 13 && HoldReasonSubCode == 2 || HoldReasonCode == 12 && HoldReasonSubCode == 2 || HoldReasonCode == 26 && HoldReasonSubCode == 120 || HoldReasonCode == 3 && HoldReasonSubCode == 0 || HoldReasonSubCode == 42)\n')
    #fsrc.write('priority = 10 \n')
    fsrc.write('+RequiresCVMFS = True \n')
    fsrc.write('Requirements = HAS_SINGULARITY == True\n')
    fsrc.write('MY.SingularityImage = "/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel9"\n')
    fsrc.write('queue '+str(n)+' from '+ifile)
    fsrc.close()

if __name__ == "__main__":
    if not len(sys.argv) > 1 or '-h' in sys.argv or '--help' in sys.argv:
        print (f"Usage: {sys.argv(0)} [-q queue] [-tree treename] [-list listfile.list] [-split S] [--sms] [--data] [--sys] [--fastsim] [--slim] [--dry-run] [--verbose] [--count] [--csv]")
        sys.exit(1)


    argv_pos    = 1
    DO_SMS      = 0
    DO_DATA     = 0
    USE_URL     = 0
    DRY_RUN     = 0
    COUNT       = 0
    VERBOSE     = 0
    CSV         = 0
    SYS         = 0
    FASTSIM     = 0
    SLIM        = 0
  
    if '-q' in sys.argv:
        p = sys.argv.index('-q')
        QUEUE = sys.argv[p+1]
        argv_pos += 2
    if '-list' in sys.argv:
        p = sys.argv.index('-list')
        LIST = sys.argv[p+1]
        argv_pos += 2
    if '-tree' in sys.argv:
        p = sys.argv.index('-tree')
        TREE = sys.argv[p+1]
        argv_pos += 2
    if '-split' in sys.argv:
        p = sys.argv.index('-split')
        SPLIT = int(sys.argv[p+1])
        argv_pos += 2
    if '--sms' in sys.argv:
        DO_SMS = 1
        argv_pos += 1
    if '--data' in sys.argv:
        DO_DATA = 1
        argv_pos += 1
    if '--dry-run' in sys.argv or '--dryrun' in sys.argv:
        DRY_RUN = 1
        argv_pos += 1
    if '--count' in sys.argv:
        COUNT = 1
        argv_pos += 1
    if '--verbose' in sys.argv:
        VERBOSE = 1
        argv_pos += 1
    if '--csv' in sys.argv:
        VERBOSE = 1
        CSV = 1
        argv_pos += 1
    if '--sys' in sys.argv:
        SYS = 1
        argv_pos += 1
    if '--fastsim' in sys.argv:
        FASTSIM = 1
        argv_pos += 1
    if '--slim' in sys.argv:
        SLIM = 1
        argv_pos += 1
        
    if SPLIT <= 1:
        SPLIT = 1
    
    print (" --- Preparing condor submission to create ntuples.")
    if DO_DATA:
        print (" --- Processing Data")

    if DO_SMS:
        print (" --- Processing SMS")
    
    if SYS:
        print (" --- Processing SYS")

    if FASTSIM:
        print (" --- Processing FastSim")

    if SLIM:
        print (" --- Processing Slim")

    if COUNT:
        print (" --- Only Counting (No Processing)")

    condor_monitor = CondorJobCountMonitor(threshold=THRESHOLD,verbose=VERBOSE)

    # input sample list
    listfile = LIST
    listname = listfile.split("/")
    listname = listname[-1]

    NAME = listname.replace(".list",'')
    
    # create and organize output folders
    TARGET  = RUN_DIR+"/"+NAME+"/"
    if not COUNT:
        os.system("rm -rf "+TARGET)
        os.system("mkdir -p "+TARGET)
    listdir = TARGET+"list/"
    srcdir  = TARGET+"src/"
    logdir  = TARGET+"log/"
    outdir  = TARGET+"out/"
    errdir  = TARGET+"err/"
    if not COUNT:
        os.system("mkdir -p "+listdir)
        os.system("mkdir -p "+logdir)
        os.system("mkdir -p "+outdir)
        os.system("mkdir -p "+errdir)
    os.system("mkdir -p "+srcdir)

    # make config directory
    config = TARGET+"config/"
    if not COUNT:
        os.system("mkdir -p "+config)

    # NOTE: there is a bug for setting the hadd verbosity level, "hadd -v 0".
    # The hadd verbosity option only works in ROOT 6.18/00 and later.
    # https://github.com/root-project/root/issues/11372
    # https://github.com/root-project/root/pull/3914

    if not COUNT:

        # make EventCount file
        if VERBOSE:
            print("making EventCount file")
        os.system("hadd "+config+"EventCount.root root/EventCount/*.root > /dev/null")
        EVTCNT = "./config/EventCount.root"

        # make FilterEff file 
        if VERBOSE:
            print("making FilterEff file")
        os.system("hadd "+config+"FilterEff.root root/FilterEff/*.root > /dev/null")
        FILTEREFF = "./config/FilterEff.root"

        # make json file
        if VERBOSE:
            print("making json file")
        os.system("cat json/GoodRunList/*.txt > "+config+"GRL_JSON.txt")
        os.system("echo -n $(tr -d '\n' < "+config+"GRL_JSON.txt) > "+config+"GRL_JSON.txt")
        JSON = "./config/GRL_JSON.txt"

        # copy xs json file
        XSJSONFILENAME = 'info_XSDB_2025-03-21_12-21.json'
        if VERBOSE:
            print("making xs json file")
        os.system(f"xrdcp -s root://cmseos.fnal.gov//store/user/z374f439/XSectionJSONs/{XSJSONFILENAME} {config}/XS_jsonfile.json")
        XSJSONFILE = f"./config/XS_jsonfile.json"

        # copy PU root files
        if VERBOSE:
            print("making Pileup file")
        os.system("cp -r root/PU "+config+".")
        PUFOLD = "./config/PU/"

        # copy BTAG SF files
        if VERBOSE:
            print("making BTAG file")
        os.system("cp -r root/BtagSF "+config+".")
        os.system("cp -r csv/BtagSF/* "+config+"BtagSF/.")
        btag_pog_fold = '/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/'
        os.system(f"cp -r {btag_pog_fold}* {config}BtagSF/")
        BTAGFOLD = "./config/BtagSF/"

        # copy LEP SF files
        if VERBOSE:
            print("making LEP file")
        os.system("cp -r root/LepSF "+config+".")
        LEPFOLD = "./config/LepSF/"

        # copy JME files
        if VERBOSE:
            print("making JME file")
        os.system("cp -r data/JME "+config+".")
        JMEFOLD = "./config/JME/"

        # copy MET trigger files
        if VERBOSE:
            print("making Trigger file")
        os.system("cp -r csv/METTrigger "+config+".")
        METFILE = "./config/METTrigger/Parameters.csv"

        # copy Prefire files
        if VERBOSE:
            print("making Prefire file")
        os.system("cp -r root/Prefire "+config+".")
        PREFIREFILE = "./config/Prefire/Prefire.root"

        if VERBOSE:
            print("Setting up working area...")
        
        os.system("cp "+EXE+" "+config+".")
        os.system("cp "+RESTFRAMES+" "+config+".")
        os.system("cp "+CMSSW_SETUP+" "+config+".")

    # output root files
    OUT_DIR = OUT_BASE+"/"+NAME+"/"
    if OUT_DIR == TARGET:
        OUT_DIR = OUT_DIR+"root/"

    total_root_files = 0

    datasetlist     = []
    clean_inputlist = []
    input_info      = {}

    # tags need to follow the format of CAMPAIGN_CMSSWX and CMSSWX must be 5 chars for later DAS checks to work
    knowntags = ["Autumn18_102X","Fall17_102X","Summer16_102X","Summer20UL16_102X","Summer20UL16APV_102X","Summer20UL17_102X","Summer20UL18_102X","Summer22_130X","Summer22EE_130X","Summer23_130X","Summer23BPix_130X"]
    
    n_samples = 0
    with open(listfile,'r') as mylist:
        inputlist = mylist.readlines()

        for flist in inputlist:
            # skip commented lines (skip if # is anywhere in line)
            if '#' in flist:
                continue

            flist = flist.strip('\n\r')

            listfile = LIST
            listname = listfile.split("/")
            listname = listname[-1]

            n_samples = n_samples+1

            dataset = flist.split("/")
            dataset = dataset[-1]
            dataset = dataset.replace(".txt",'')

            filetag = ""
            for ktag in knowntags:
                if ktag in flist:
                    filetag = ktag

            # get list of ROOT files 
            rootlist = []
            with open(flist,'r') as myflist:
                inputfilelist = myflist.readlines()

                for afile in inputfilelist:
                    afile = afile.strip('\n\r')
                    rootlist.append(afile);
            
            n_root_files        = len(rootlist)
            n_jobs              = SPLIT * n_root_files
            total_root_files    += n_root_files
            
            input_info[dataset+'_'+filetag] = {}
            clean_inputlist.append(dataset+'_'+filetag)
            input_info[dataset+'_'+filetag]["n_root_files"]   = n_root_files
            input_info[dataset+'_'+filetag]["n_jobs"]         = n_jobs

            if len(datasetlist) == 0:
                datasetlist.append((dataset,filetag,rootlist))
                if not COUNT:
                    os.system("rm -rf "+OUT_DIR+dataset+"_"+filetag+"/")
                    os.system("mkdir -p "+OUT_DIR+dataset+"_"+filetag+"/")
                continue
            
            tagtuple = [item for item in datasetlist if item[0] == dataset]
            if len(tagtuple) == 0:
                datasetlist.append((dataset,filetag,rootlist))
                if not COUNT:
                    os.system("rm -rf "+OUT_DIR+dataset+"_"+filetag+"/")
                    os.system("mkdir -p "+OUT_DIR+dataset+"_"+filetag+"/")
                continue

            p = datasetlist.index(tagtuple[0])
            datasetlist[p][2].extend(rootlist)

    if VERBOSE:
        print("Created area for output root files")

    for (dataset,filetag,rootlist) in datasetlist:
        if not COUNT:
            os.system("mkdir -p "+os.path.join(listdir, dataset+'_'+filetag))
            listdir_sam = os.path.join(listdir, dataset+'_'+filetag)
            listlist = create_filelist(rootlist, dataset, filetag)
            overlist_name = listdir_sam+'/'+dataset+'_'+filetag+'_list.list'
            with open(overlist_name,'w') as overlist:
                newlistlist = ['config/'+'/'.join(l.split('/')[-3:])+'\n' for l in listlist]
                overlist.writelines(newlistlist)
                overlist.close()

        if not COUNT:
            os.system("mkdir -p "+os.path.join(logdir, dataset+'_'+filetag))
            os.system("mkdir -p "+os.path.join(outdir, dataset+'_'+filetag))
            os.system("mkdir -p "+os.path.join(errdir, dataset+'_'+filetag))

            file_name = os.path.join(OUT_DIR, dataset+'_'+filetag, overlist_name.split('/')[-1].replace('_list.list', '_$(ItemIndex)_$(Step)'))

            logfile = os.path.join(logdir, dataset+'_'+filetag, file_name.split('/')[-1])
            outfile = os.path.join(outdir, dataset+'_'+filetag, file_name.split('/')[-1])
            errfile = os.path.join(errdir, dataset+'_'+filetag, file_name.split('/')[-1])

            script_name = srcdir+'_'.join([dataset, filetag])+'.submit'
            write_sh(script_name, overlist_name, file_name+'.root', logfile, outfile, errfile, dataset, filetag, SPLIT, NAME)
            write_sh_single(script_name, overlist_name, file_name+'.root', logfile, outfile, errfile, dataset, filetag, SPLIT, NAME)
    
    if VERBOSE:
        print("Created area for log files")

    if not COUNT:
        #print listdir
        os.system("cp -r "+listdir+" "+config)
        #print "creating tarball from: ", TARGET
        os.system("sleep 10") # sleep so copy command(s) can catch up...
        os.system("tar -C "+config+"/../ -czf "+TARGET+"/config.tgz config")
        # Why was this needed...???
        # os.system("mkdir -p /ospool/cms-user/"+USER+"/"+NAME)
        # os.system("cp "+TARGET+"/config.tgz /ospool/cms-user/"+USER+"/"+NAME+"/config.tgz")
        if VERBOSE:
            print("Created tar ball")

    submit_dir  = srcdir        
    submit_list = [os.path.join(submit_dir, f) for f in os.listdir(submit_dir) if (os.path.isfile(os.path.join(submit_dir, f)) and ('.submit' in f) and ('_single' not in f))]
    #n_samples   = len(submit_list)
    total_jobs  = SPLIT * total_root_files

    # Prep csv file
    if CSV:
        csv_name = LIST.split("/")[-1].split(".")[0]
        f_csv = open(TARGET+"/"+csv_name+".csv",'w')
        f_csv.write('sample,clusterid,totaljobs')
        f_csv.write('\n')

    # don't submit jobs if --dry-run or --count are used
    if not DRY_RUN and not COUNT:
        for f in submit_list:
            sample_handle = f.split("/")
            sample_handle = sample_handle[-1]
            sample_handle = sample_handle.replace(".submit",'')
            print (f"submitting: {f}")
            if CSV:
                condor_monitor.wait_until_jobs_below()
                os.system('condor_submit '+f+' | tee '+sample_handle+'.txt')
                with open(sample_handle+'.txt','r') as sample_submit_file:
                    lines = sample_submit_file.read()
                    clusterid_index = lines.find('cluster ')
                    if clusterid_index != -1:
                        clusterid = lines[clusterid_index+len('cluster '):]
                        input_info[sample_handle]["clusterid"] = clusterid
                os.system('rm '+sample_handle+'.txt')
            else:
                condor_monitor.wait_until_jobs_below()
                os.system('condor_submit ' + f)
            
    
    # Number of ROOT files and jobs per sample 
    if VERBOSE:
        for f in clean_inputlist:
            n_root_files    = input_info[f]["n_root_files"] 
            n_jobs          = input_info[f]["n_jobs"] 
            n_jobs = SPLIT * n_root_files
            print(f"sample: {f}")
            print(f" - number of root files  = {n_root_files}")
            print(f" - number of jobs        = {n_jobs}")
            # make sure that "clusterid" has been filled to avoid key error
            if not DRY_RUN and not COUNT and CSV:
                f_csv.write(f"{0}".format(f+","+input_info[f]["clusterid"].replace('.\n','')+",{0}".format(n_jobs)+'\n'))

    # Close csv file
    if CSV:
         f_csv.close()
    
    # Summary Info
    print ("----------------------------")
    print ("Condor Submission Info")
    print ("----------------------------")
    print (f"sample list:             {LIST}")
    print (f"working directory:       {TARGET}")
    print (f"output directory:        {OUT_DIR}")
    print (f"number of samples:       {n_samples}")
    print (f"split:                   {SPLIT}")
    print (f"total input root files:  {total_root_files}")
    print (f"total condor jobs:       {total_jobs}")
    print ("----------------------------")

    if DRY_RUN or COUNT:
        print("No jobs were submitted.")
    else:
        print(Fore.GREEN + "Congrats... your jobs were submitted!" + Fore.RESET)

