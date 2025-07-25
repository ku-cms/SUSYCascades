#! /usr/bin/env python
import os, sys, time, subprocess
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'python')))
from CondorJobCountMonitor import CondorJobCountMonitor

# ----------------------------------------------------------- #
# Parameters
# ----------------------------------------------------------- #
# current working directory
pwd         = os.environ['PWD']
RUN_DIR     = pwd
TEMP        = pwd
jobEXE      = "execute_script_EventCount.sh"
EXE         = "MakeEventCount_NANO.x"
RESTFRAMES  = './scripts/setup_RestFrames_connect.sh'
#CMSSW_SETUP = './scripts/cmssw_setup_connect.sh'
CMSSW_SETUP = './scripts/cmssw_setup_connect_el9.sh'
TREE        = "Runs"
USER        = os.environ['USER']
OUT         = "/uscms/home/"+USER+"/nobackup/EventCount/root/"
#OUT         = "/ospool/cms-user/"+USER+"/NTUPLES/Processing"
LIST        = "default.list"
QUEUE       = ""
MAXN        = 1
CONNECT     = False
THRESHOLD   = 90000
# ----------------------------------------------------------- #

def get_auto_THRESHOLD():
    result = subprocess.run(
        ["condor_config_val", "-dump"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )
    for line in result.stdout.splitlines():
        if line.startswith("MAX_JOBS_PER_OWNER"):
            _, value = line.split("=")
            return int(value.strip())
            break
    return MAX_JOBS_SUB

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
        if len(sublist) >= MAXN and MAXN > 0:
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

def write_sh(srcfile,ifile,ofile,logfile,outfile,errfile,dataset,filetag):
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
    if DO_CASCADES == 1:
        fsrc.write('--cascades ')
    if DO_PRIVATEMC == 1:
        fsrc.write('--privateMC ')
    fsrc.write('-dataset='+dataset+" ")
    fsrc.write('-filetag='+filetag+" ")

    fsrc.write('\n')
    outlog = outfile+".out"
    errlog = errfile+".err"
    loglog = logfile+".log"
    fsrc.write('output = '+outlog+" \n")
    fsrc.write('error = '+errlog+" \n")
    fsrc.write('log = '+loglog+" \n")
    fsrc.write('priority = 10 \n')
    fsrc.write('Requirements = (Machine != "red-node000.unl.edu" && Machine != "ncm*.hpc.itc.rwth-aachen.de" && Machine != "mh-epyc7662-8.t2.ucsd.edu")\n')
    fsrc.write('request_memory = 1 GB \n')

    if CONNECT is True:
        transfer_input = 'transfer_input_files = '+TARGET+'config.tgz,/ospool/cms-user/zflowers/public/sandbox-CMSSW_13_3_1-el9.tar.bz2\n'
    else:
        transfer_input = 'transfer_input_files = '+TARGET+'config.tgz,/uscms/home/z374f439/nobackup/whatever_you_want/sandbox-CMSSW_13_3_1-el9.tar.bz2\n'
    fsrc.write(transfer_input)

    fsrc.write('should_transfer_files = YES\n')
    fsrc.write('when_to_transfer_output = ON_EXIT\n')

    transfer_out_files = 'transfer_output_files = '+ofile.split('/')[-1]+'\n'
    fsrc.write(transfer_out_files)

    transfer_out_remap = 'transfer_output_remaps = "'+ofile.split('/')[-1]+'='+ofile+'"\n'
    fsrc.write(transfer_out_remap)
    
    fsrc.write('+ProjectName="cms.org.ku"\n')
    fsrc.write('+RequiresCVMFS = True \n')
    fsrc.write('Requirements = HAS_SINGULARITY == True\n')
    fsrc.write('RequestCpus = 1\n')
    fsrc.write('periodic_hold_subcode = 42\n')
    fsrc.write('periodic_hold = (CpusUsage > RequestCpus) && (JobStatus == 2) && (CurrentTime - EnteredCurrentStatus > 60)\n')
    fsrc.write('+JobTransforms = "if HoldReasonSubCode == 42 set RequestCpus = MIN(RequestCpus + 1, 32)"\n')
    fsrc.write('periodic_release = (HoldReasonCode == 12 && HoldReasonSubCode == 256 || HoldReasonCode == 13 && HoldReasonSubCode == 2 || HoldReasonCode == 12 && HoldReasonSubCode == 2 || HoldReasonCode == 26 && HoldReasonSubCode == 120 || HoldReasonCode == 3 && HoldReasonSubCode == 0 || HoldReasonSubCode == 42)\n')
    fsrc.write('+REQUIRED_OS="rhel9"\n')
    fsrc.write('job_lease_duration = 3600\n')
    fsrc.write('on_exit_remove = (ExitCode == 0) || (NumJobStarts >= 3)\n')
    fsrc.write('max_retries = 3\n')
    fsrc.write('MY.SingularityImage = "/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel9"\n')
    fsrc.write('queue from '+ifile+'\n')
    fsrc.close()

def write_sh_single(srcfile,ifile,ofile,logfile,outfile,errfile,dataset,filetag):
    srcfile = srcfile.replace('.submit','_single.submit')
    ofile = ofile.replace('_$(ItemIndex)','_0')
    outfile = outfile.replace('_$(ItemIndex)','_0')
    errfile = errfile.replace('_$(ItemIndex)','_0')
    logfile = logfile.replace('_$(ItemIndex)','_0')
    list_name = os.path.basename(ifile).replace('_list', '_0')
    ifile = os.path.join('./config/list', list_name)

    fsrc = open(srcfile,'w')
    fsrc.write('# Note: For only submitting 1 job! \n')
    fsrc.write('# output and list args need to be updated \n')
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
    if DO_CASCADES == 1:
        fsrc.write('--cascades ')
    if DO_PRIVATEMC == 1:
        fsrc.write('--privateMC ')
    fsrc.write('-dataset='+dataset+" ")
    fsrc.write('-filetag='+filetag+" ")
    fsrc.write('\n')
    outlog = outfile+".out"
    errlog = errfile+".err"
    loglog = logfile+".log"
    fsrc.write('output = '+outlog+" \n")
    fsrc.write('error = '+errlog+" \n")
    fsrc.write('log = '+loglog+" \n")
    fsrc.write('request_memory = 2 GB \n')
    fsrc.write('+RequiresCVMFS = True \n')

    if CONNECT is True:
        transfer_input = 'transfer_input_files = '+TARGET+'config.tgz,/ospool/cms-user/zflowers/public/sandbox-CMSSW_13_3_1-el9.tar.bz2\n'
    else:
        transfer_input = 'transfer_input_files = '+TARGET+'config.tgz,/uscms/home/z374f439/nobackup/whatever_you_want/sandbox-CMSSW_13_3_1-el9.tar.bz2\n'
    fsrc.write(transfer_input)

    fsrc.write('should_transfer_files = YES\n')
    fsrc.write('when_to_transfer_output = ON_EXIT\n')

    transfer_out_files = 'transfer_output_files = '+ofile.split('/')[-1]+'\n'
    fsrc.write(transfer_out_files)

    transfer_out_remap = 'transfer_output_remaps = "'+ofile.split('/')[-1]+'='+ofile+'"\n'
    fsrc.write(transfer_out_remap)
    
    fsrc.write('+ProjectName="cms.org.ku"\n')
    fsrc.write('+RequiresCVMFS = True \n')
    fsrc.write('Requirements = HAS_SINGULARITY == True\n')
    fsrc.write('RequestCpus = 1\n')
    fsrc.write('periodic_hold_subcode = 42\n')
    fsrc.write('periodic_hold = (CpusUsage > RequestCpus) && (JobStatus == 2) && (CurrentTime - EnteredCurrentStatus > 60)\n')
    fsrc.write('+JobTransforms = "if HoldReasonSubCode == 42 set RequestCpus = MIN(RequestCpus + 1, 32)"\n')
    fsrc.write('periodic_release = (HoldReasonCode == 12 && HoldReasonSubCode == 256 || HoldReasonCode == 13 && HoldReasonSubCode == 2 || HoldReasonCode == 12 && HoldReasonSubCode == 2 || HoldReasonCode == 26 && HoldReasonSubCode == 120 || HoldReasonCode == 3 && HoldReasonSubCode == 0 || HoldReasonSubCode == 42)\n')
    fsrc.write('+REQUIRED_OS="rhel9"\n')
    fsrc.write('job_lease_duration = 3600\n')
    fsrc.write('on_exit_remove = (ExitCode == 0) || (NumJobStarts >= 3)\n')
    fsrc.write('max_retries = 3\n')
    fsrc.write('MY.SingularityImage = "/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel9"\n')
    fsrc.write('queue\n')
    fsrc.close()

if __name__ == "__main__":
    if not len(sys.argv) > 1 or '-h' in sys.argv or '--help' in sys.argv:
        print (f"Usage: {sys.argv(0)} [-q queue] [-tree treename] [-list listfile.list] [--sms] [--connect]")
        sys.exit(1)

    argv_pos     = 1
    DO_SMS       = 0
    DO_DATA      = 0
    DRY_RUN      = 0
    DO_CASCADES  = 0
    DO_PRIVATEMC = 0
  
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
    if '--sms' in sys.argv:
        DO_SMS = 1
        argv_pos += 1
        TREE = "Events"
    if '--cascades' in sys.argv:
        DO_CASCADES = 1
        argv_pos += 1
    if '--private' in sys.argv:
        DO_PRIVATEMC = 1
        argv_pos += 1
    if '--privateMC' in sys.argv:
        DO_PRIVATEMC = 1
        argv_pos += 1
    if '--connect' in sys.argv:
        CONNECT = True
        argv_pos += 1
    if '--dry-run' in sys.argv or '--dryrun' in sys.argv:
        DRY_RUN = 1
        argv_pos += 1
    
    if DO_DATA:
        print ("Processing Data")

    if DO_SMS:
        print ("Processing as SMS")

    if DO_CASCADES:
        print ("Processing as Cascades")

    if DO_PRIVATEMC:
        print ("Processing as Private MC")

    THRESHOLD = get_auto_THRESHOLD()

    # input sample list
    listfile = LIST
    listname = listfile.split("/")
    listname = listname[-1]

    NAME = listname.replace(".list",'')
    NAME += "_EventCount"
    
    print(listname)
    print(NAME)
    print(RUN_DIR)
        
    # create and organize output folders
    TARGET  = RUN_DIR+"/"+NAME+"/"
    os.system("rm -rf "+TARGET)
    os.system("mkdir -p "+TARGET)
    listdir = TARGET+"list/"
    srcdir  = TARGET+"src/"
    logdir  = TARGET+"log/"
    outdir  = TARGET+"out/"
    errdir  = TARGET+"err/"
    os.system("mkdir -p "+listdir)
    os.system("mkdir -p "+logdir)
    os.system("mkdir -p "+outdir)
    os.system("mkdir -p "+errdir)
    os.system("mkdir -p "+srcdir)

    # make config directory
    config = TARGET+"config/"
    os.system("mkdir -p "+config)

    os.system("cp "+EXE+" "+config+".")
    os.system("cp "+RESTFRAMES+" "+config+".")
    os.system("cp "+CMSSW_SETUP+" "+config+".")

    if CONNECT is True:
        OUT  = "/ospool/cms-user/"+USER+"/EventCount/root/"

    # output root files
    ROOT = OUT+"/"+NAME+"/"
    if ROOT == TARGET:
        ROOT = ROOT+"root/"

    datasetlist = []

    # tags need to follow the format of CAMPAIGN_CMSSWX and CMSSWX must be 5 chars for later DAS checks to work
    knowntags = ["Autumn18_102X","Fall17_102X","Summer16_102X","Summer20UL16_106X","Summer20UL16APV_106X","Summer20UL17_106X","Summer20UL18_106X","Summer22_130X","Summer22EE_130X","Summer23_130X","Summer23BPix_130X"]
    
    with open(listfile,'r') as mylist:
        inputlist = mylist.readlines()

        for flist in inputlist:
            # skip commented lines (skip if # is anywhere in line)
            if '#' in flist:
                continue
            flist = flist.strip('\n\r')
            print(f"Processing list from {flist}")

            listfile = LIST
            listname = listfile.split("/")
            listname = listname[-1]

            dataset = flist.split("/")
            dataset = dataset[-1]
            dataset = dataset.replace(".txt",'')

            # make output folders

            filetag = ""
            for ktag in knowntags:
                if ktag in flist:
                    filetag = ktag

            rootlist = []
            with open(flist,'r') as myflist:
                inputfilelist = myflist.readlines();

                for afile in inputfilelist:
                    afile = afile.strip('\n\r')
                    rootlist.append(afile);

            if len(datasetlist) == 0:
                datasetlist.append((dataset,filetag,rootlist))
                os.system("rm -rf "+ROOT+dataset+"_"+filetag+"/")
                os.system("mkdir -p "+ROOT+dataset+"_"+filetag+"/")
                continue
            
            tagtuple = [item for item in datasetlist if item[0] == dataset]
            if len(tagtuple) == 0:
                datasetlist.append((dataset,filetag,rootlist))
                os.system("mkdir -p "+ROOT+dataset+"_"+filetag+"/")
                continue

            p = datasetlist.index(tagtuple[0])
            datasetlist[p][2].extend(rootlist)

    job_total_dataset = {}
    for (dataset,filetag,rootlist) in datasetlist:
        os.system("mkdir -p "+os.path.join(listdir, dataset+'_'+filetag))
        listdir_sam = os.path.join(listdir, dataset+'_'+filetag)
        listlist = create_filelist(rootlist, dataset, filetag)
        overlist_name = listdir_sam+'/'+dataset+'_'+filetag+'_list.list'
        with open(overlist_name,'w') as overlist:
            newlistlist = ['config/'+'/'.join(l.split('/')[-3:])+'\n' for l in listlist]
            overlist.writelines(newlistlist)
            overlist.close()

        os.system("mkdir -p "+os.path.join(logdir, dataset+'_'+filetag))
        os.system("mkdir -p "+os.path.join(outdir, dataset+'_'+filetag))
        os.system("mkdir -p "+os.path.join(errdir, dataset+'_'+filetag))

        file_name = os.path.join(ROOT, dataset+'_'+filetag, overlist_name.split('/')[-1].replace('_list.list', '_$(ItemIndex)'))

        logfile = os.path.join(logdir, dataset+'_'+filetag, file_name.split('/')[-1])
        outfile = os.path.join(outdir, dataset+'_'+filetag, file_name.split('/')[-1])
        errfile = os.path.join(errdir, dataset+'_'+filetag, file_name.split('/')[-1])

        script_name = srcdir+'_'.join([dataset, filetag])+'.submit'
        write_sh(script_name, overlist_name, file_name+'.root', logfile, outfile, errfile, dataset, filetag)
        write_sh_single(script_name, overlist_name, file_name+'.root', logfile, outfile, errfile, dataset, filetag)
        job_total_dataset[dataset+'_'+filetag] = len(rootlist)

    #print listdir
    os.system("cp -r "+listdir+" "+config)
    os.system("tar -C "+config+"/../ -czf "+TARGET+"/config.tgz config")

    submit_dir  = srcdir        
    submit_list = [os.path.join(submit_dir, f) for f in os.listdir(submit_dir) if (os.path.isfile(os.path.join(submit_dir, f)) and ('.submit' in f) and ('_single' not in f))]

    if not DRY_RUN:
        for f in submit_list:
            sample_handle = f.split("/")
            sample_handle = sample_handle[-1]
            sample_handle = sample_handle.replace(".submit",'')
            print (f"submitting: {f}")
            submit_jobs = job_total_dataset[sample_handle]
            condor_monitor = CondorJobCountMonitor(threshold=THRESHOLD-submit_jobs,verbose=False)
            condor_monitor.wait_until_jobs_below()
            print("submitting: ", f)
            os.system('condor_submit ' + f)
   
