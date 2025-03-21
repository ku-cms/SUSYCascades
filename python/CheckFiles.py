# ------------------------------------------------------------------------------------------------
# Adopted from https://github.com/ku-cms/SUSYCascades/blob/master/MakeNtuple_Connect/CheckJobs.py
# Example call of script from SUSYCascades:
#    python3 python/CheckFiles.py -d Summer23_130X/ -o ../../../NTUPLES/Processing/Summer23_130X/ -e -r
# ------------------------------------------------------------------------------------------------
import os, argparse, subprocess, itertools, ROOT
import time
from new_countEvents import EventCount as EventCount
from CondorJobCountMonitor import CondorJobCountMonitor
USER = os.environ['USER']
DO_EVENTCOUNT = False # For checking event count files

def getMissingFilesEC(outputDir,nList):
    baseIndexList = list(range(nList))
    for filename in os.listdir(outputDir):
        if not filename.endswith('root'): continue;
        filename_str = filename.split(".")[0]
        file_index = int(filename_str.split("_")[-1])
        if file_index in baseIndexList:
            baseIndexList.remove(file_index)
    return baseIndexList

def getMissingFiles(outputDir,nSplit,nList):
    if DO_EVENTCOUNT:
        return getMissingFilesEC(outputDir,nList)
    baseTupleList = list(itertools.product(range(nList), range(nSplit)))
    for filename in os.listdir(outputDir):
        if not filename.endswith('root'): continue;
        filename_str = filename.split(".")[0]
        Tuple = (int(filename_str.split("_")[-2]),int(filename_str.split("_")[-1]))
        if Tuple in baseTupleList:
            baseTupleList.remove(Tuple)
    return baseTupleList

def makeSubmitScript(tuple_pairs,submitName,resubmit,maxResub,DataSetName):
    tuple_filelist = f"{submitName}_tuple.txt"
    resubmitFiles = 0
    if tuple_pairs == []:
        return # Have no jobs to resubmit
    tuple_pairs = list(set(tuple_pairs)) # remove potential dupes
    with open(tuple_filelist,'w') as tuple_file:
        for tuple_pair in tuple_pairs:
            resubmitFiles = resubmitFiles+1
            if DO_EVENTCOUNT: tuple_file.write(f"{tuple_pair}\n")
            else: tuple_file.write(f"{tuple_pair[0]},{tuple_pair[1]}\n")
    newFileName = f"{submitName}_tuple.submit"
    os.system(f"cp {submitName}_single.submit {newFileName}")
    with open(newFileName,'r') as file:
        file_content = file.read()
    if DO_EVENTCOUNT:
        file_content = file_content.replace("0.list","$(list).list")
        file_content = file_content.replace("0.root","$(list).root")
        file_content = file_content.replace("0.out","$(list).out")
        file_content = file_content.replace("0.err","$(list).err")
        file_content = file_content.replace("0.log","$(list).log")
        file_content = file_content.replace("queue",f"queue list from {tuple_filelist}")
    else:
        file_content = file_content.replace("0_0","$(list)_$(split)")
        file_content = file_content.replace("-split=1","-split=$$([$(split)+1])")
        file_content = file_content.replace("X_0.list","X_$(list).list")
        file_content = file_content.replace("queue",f"queue list,split from {tuple_filelist}")
    file_content = file_content.replace('RequestCpus = 1','RequestCpus = 4')
    with open(newFileName, 'w') as file:
        file.write(file_content)
    # Check number of resubmitFiles
    print('total resubmit:',resubmitFiles)
    if resubmit:
        if resubmitFiles > maxResub:
            print(f"You are about to make {resubmitFiles} and resubmit {resubmitFiles} jobs for dataset: {DataSetName}!")
            print(f"You should double check there are no issues with your condor submissions.")
            print(f"If you are confident you want to resubmit, then you should rerun this script with '-l {resubmitFiles}'.")
        else:
            condor_monitor = CondorJobCountMonitor(threshold=80000,verbose=False)
            condor_monitor.wait_until_jobs_below()
            os.system(f"condor_submit {newFileName}")

def testRootFile(root_file):
    status = True
    try:
        # Can ignore any message like: ReadStreamerInfo, class:string, illegal uid=-2
        # https://root-forum.cern.ch/t/readstreamerinfo-illegal-uid-with-newer-root-versions/41073
        ROOT.gErrorIgnoreLevel = ROOT.kFatal+1
        root_file_test = ROOT.TFile.Open(root_file);
        if not root_file_test or root_file_test.IsZombie():
            status = False
        if root_file_test.IsOpen():
            root_file_test.Close()
    except:
        status = False
    return status

def eventCountCheck(directory):
    failed_nums = []
    event_count = EventCount()
    for dirpath, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.root'):
                if not event_count.checkInternalDASCount(os.path.join(dirpath, file)):
                    file = file.split(".root")[0]
                    if DO_EVENTCOUNT: failed_nums.append(int(file.split("_")[-1]))
                    else: failed_nums.append((int(file.split("_")[-2]), int(file.split("_")[-1])))
    return failed_nums

# Check condor jobs
def checkJobs(workingDir,outputDir,skipEC,skipDAS,skipMissing,skipSmall,skipErr,skipOut,skipZombie,resubmit,maxResub,filter_list):
    print("Running over the directory '{0}'.".format(workingDir))
    print("------------------------------------------------------------")
    grep_ignore = "-e \"Warning\" -e \"WARNING\" -e \"TTree::SetBranchStatus\" -e \"libXrdSecztn.so\" -e \"Phi_mpi_pi\" -e \"tar: stdout: write error\" -e \"INFO\""
    grep_ignore += " -e \"SetBranchAddress\""
    srcDir = os.listdir(workingDir+"/src/")
    nJobs = 0
    for file in srcDir:
        if "X.submit" not in file:
            continue
        do_name = False
        if filter_list:
            for name in filter_list:
                if name in file:
                    do_name = True
                    break
        if filter_list and not do_name:
            continue
        DataSetName = file.split(".")[0]
        resubmitFiles = []
        if(not skipMissing):
            with open(workingDir+"/src/"+DataSetName+".submit") as file:
                file.seek(0,2)
                file.seek(file.tell()-2,0)
                while file.read(1) != "\n" and file.tell() > 0:
                    file.seek(file.tell()-2,0)
                last_line = file.readline().strip()
                nSplit = 1
                if not DO_EVENTCOUNT:
                    nSplit = int(last_line.split(" ")[1])
            listDir = os.path.abspath(workingDir+"/list/"+DataSetName+"/")
            listFiles = os.listdir(listDir)
            numList = len([listFile for listFile in listFiles if os.path.isfile(os.path.join(listDir,listFile))])
            numList = numList - 1 # remove master list file
            nJobsSubmit = numList*nSplit
            resubmitFiles.extend(getMissingFiles(outputDir+'/'+DataSetName,nSplit,numList))
            print(f"Got {len(resubmitFiles)} missing files for dataset {DataSetName}")
        if(not skipSmall):
            # define size cutoff using -size: find will get all files below the cutoff
            bash = "find "+outputDir+'/'+DataSetName+" -type f -size -1k"
            smallFiles = subprocess.check_output(['bash','-c', bash]).decode()
            smallFiles = smallFiles.split("\n")
            smallFiles.remove('')
            num_small = 0
            for smallFile in smallFiles:
                smallFile = smallFile.split(".root")[0]
                Tuple = None
                if DO_EVENTCOUNT: Tuple = int(smallFile.split("_")[-1])
                else: Tuple = (int(smallFile.split("_")[-2]),int(smallFile.split("_")[-1]))
                resubmitFiles.append(Tuple)
                num_small = num_small+1
            print(f"Got {num_small} small files for dataset {DataSetName}")
        if(not skipErr):
            bash = "grep -r -v "+grep_ignore+" "+workingDir+"/err/"+DataSetName+"/"
            # "errorFiles" is actually all lines with errors that we don't ignore
            # thus, the number of "errorFiles" is often much greater than number of files with an actual error
            errorFiles = subprocess.check_output(['bash','-c',bash]).decode()
            errorFiles = errorFiles.split("\n")
            errorFiles.remove('')
            num_error = 0
            # loop over all lines with errors
            for errorFile in errorFiles:
                errorFile = errorFile.split(".err")[0]
                Tuple = None
                if DO_EVENTCOUNT: Tuple = int(errorFile.split("_")[-1])
                else: Tuple = (int(errorFile.split("_")[-2]),int(errorFile.split("_")[-1]))
                if Tuple not in resubmitFiles:
                    # count number of files with an error
                    num_error = num_error+1
                    resubmitFiles.append(Tuple)
            print(f"Got {num_error} error files for dataset",DataSetName)
        if(not skipOut):
            outFiles = []
            bash = "grep -r \"Ntree 0\" "+ workingDir +"/out/"+DataSetName+"/"
            try:
                outFiles = subprocess.check_output(['bash','-c',bash]).decode()
                outFiles = outFiles.split("\n")
                outFiles.remove('')
            except subprocess.CalledProcessError as e:
                pass # exception means no out files have issues so just keep going
            num_out = 0
            for outFile in outFiles:
                num_out = num_out+1
                outFile = outFile.split(".out")[0]
                Tuple = None
                if DO_EVENTCOUNT: Tuple = int(outFile.split("_")[-1])
                else: Tuple = (int(outFile.split("_")[-2]),int(outFile.split("_")[-1]))
                if Tuple not in resubmitFiles:
                    resubmitFiles.append(Tuple)
            print(f"Got {num_out} out files for dataset",DataSetName)
        if(not skipZombie):
            num_zomb = 0
            bash = "find "+outputDir+'/'+DataSetName+" -type f"
            zombFiles = subprocess.check_output(['bash','-c', bash]).decode()
            zombFiles = zombFiles.split("\n")
            zombFiles.remove('')
            for zombFile in zombFiles:
                if testRootFile(zombFile) is True:
                    continue
                zombFile = zombFile.split(".root")[0]
                Tuple = None
                if DO_EVENTCOUNT: Tuple = int(zombFile.split("_")[-1])
                else: Tuple = (int(zombFile.split("_")[-2]),int(zombFile.split("_")[-1]))
                resubmitFiles.append(Tuple)
                num_zomb = num_zomb+1
            print(f"Got {num_zomb} zombie root files for dataset",DataSetName)
        if(not skipEC):
            failedEC = eventCountCheck(outputDir+'/'+DataSetName)
            resubmitFiles.extend(failedEC)
            print(f"Got {len(failedEC)} files failed EventCount check for dataset",DataSetName)
        if(not skipDAS):
            NDAS_true = 0
            NDAS = 0
            event_count = EventCount()
            listDir = os.path.join(workingDir, "config/list", DataSetName)
            for lFile in os.listdir(listDir):
                if lFile.endswith('.list') and not lFile.endswith('_list.list'):
                    with open(os.path.join(listDir,lFile),'r') as input_root_files:
                        input_root_filename = input_root_files.readline().strip()
                        dataset = event_count.GetDatasetFromFile(input_root_filename)
                        NDAS_true = event_count.EventsInDAS(dataset, False)
                        break
                        # Can get total per file if needed
                        #NDAS_true += event_count.EventsInDAS(input_root_filename) # get true total events from DAS
            outfiles = os.listdir(outputDir+'/'+DataSetName)
            for outfile in outfiles:
                NDAS += event_count.GetDASCount(os.path.join(outputDir+'/'+DataSetName,outfile))
            if NDAS != NDAS_true:
                print(f'{DataSetName} failed the DAS check! ({100.*round(NDAS/NDAS_true,3)}%) Use other options to investigate')
            else:
                print(f'{DataSetName} passed the DAS check!')
        makeSubmitScript(resubmitFiles,workingDir+"/src/"+DataSetName,resubmit,maxResub,DataSetName)
        nJobs += len(resubmitFiles)
    return nJobs

def main():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--eventCount",           action='store_true', help="for checking event count jobs")
    parser.add_argument("--directory",      "-d", default="Summer22_130X/", help="working directory for condor submission")
    parser.add_argument("--output",         "-o", default="/ospool/cms-user/"+USER+"/NTUPLES/Processing/Summer22_130X/", help="output area for root files")
    parser.add_argument("--skipResubmit",   "-r", action='store_false', help="do not resubmit files")
    parser.add_argument("--skipEC",         "-c", action='store_true', help="skip checking event count tree")
    parser.add_argument("--skipDAS",        "-w", action='store_true', help="skip checking compared to central DAS")
    parser.add_argument("--skipMissing",    "-m", action='store_true', help="skip checking missing files")
    parser.add_argument("--skipSmall",      "-s", action='store_true', help="skip checking small files")
    parser.add_argument("--skipErr",        "-e", action='store_true', help="skip checking err files (recommended: slow)")
    parser.add_argument("--skipOut",        "-u", action='store_true', help="skip checking out files")
    parser.add_argument("--skipZombie",     "-z", action='store_true', help="skip checking zombie root files")
    parser.add_argument("--maxResub",       "-l", default=100, help="max number of jobs to resubmit")

    global DO_EVENTCOUNT
    options        = parser.parse_args()
    DO_EVENTCOUNT  = options.eventCount
    directory      = options.directory	
    output         = options.output
    resubmit       = options.skipResubmit
    skipEC         = options.skipEC
    skipDAS        = options.skipDAS
    skipMissing    = options.skipMissing
    skipSmall      = options.skipSmall
    skipErr        = options.skipErr
    skipOut        = options.skipOut
    skipZombie     = options.skipZombie
    maxResub       = int(options.maxResub)

    filter_list = [ # list of name of samples to explicitly check (empty list will check all)
        #"WZ_TuneCP5_13p6TeV_pythia8",
    ]

    nJobs = checkJobs(directory,output,skipEC,skipDAS,skipMissing,skipSmall,skipErr,skipOut,skipZombie,resubmit,maxResub,filter_list)

if __name__ == "__main__":
    main()
