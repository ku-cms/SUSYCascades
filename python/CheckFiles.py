# ------------------------------------------------------------------------------------------------
# Adopted from https://github.com/ku-cms/SUSYCascades/blob/master/MakeNtuple_Connect/CheckJobs.py
# Example call of script from SUSYCascades:
#    python3 python/CheckFiles.py -d Summer23_130X/ -o ../../../NTUPLES/Processing/Summer23_130X/ -e -r
# ------------------------------------------------------------------------------------------------
import os, argparse, subprocess, itertools, ROOT, time, math, re, glob
from collections import defaultdict
from new_countEvents import EventCount as EventCount
from CondorJobCountMonitor import CondorJobCountMonitor
USER = os.environ['USER']
DO_EVENTCOUNT = False # For checking event count files

def has_rule(did, rse):
    cmd = ["rucio", "rule", "list", "--did", did]
    try:
        out = subprocess.check_output(cmd, text=True)
    except subprocess.CalledProcessError:
        return False
    return rse in out

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

def getMissingFilesEC(outputDir,nList):
    baseIndexList = list(range(nList))
    for filename in os.listdir(outputDir):
        if not filename.endswith('root'): continue;
        filename_str = filename.split(".")[0]
        file_index = int(filename_str.split("_")[-1])
        if file_index in baseIndexList:
            baseIndexList.remove(file_index)
    return baseIndexList

def getMissingFiles(outputDir, listFile):
    """
    Returns a sorted list of missing (fileIndex, splitIndex) tuples
    using the master list as truth and checking the outputDir.
    """
    # Build a map: fileIndex -> set of expected splitIndices
    expectedMap = defaultdict(set)
    with open(listFile) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            try:
                fileIndex = int(parts[1])
                splitIndex = int(parts[2])
                expectedMap[fileIndex].add(splitIndex)
            except ValueError:
                continue

    # Remove found splitIndices from the map
    for filename in os.listdir(outputDir):
        if not filename.endswith('.root'):
            continue
        parts = filename.split(".")[0].split("_")
        try:
            fileIndex = int(parts[-2])
            splitIndex = int(parts[-1])
            if fileIndex in expectedMap:
                expectedMap[fileIndex].discard(splitIndex)
                if not expectedMap[fileIndex]:
                    # Remove fileIndex if all splits are accounted for
                    del expectedMap[fileIndex]
        except (IndexError, ValueError):
            continue

    # Flatten the map into a sorted list of tuples
    missingTuples = []
    for fileIndex, splits in expectedMap.items():
        for splitIndex in splits:
            missingTuples.append((fileIndex, splitIndex))

    return sorted(missingTuples)

# helper to find the master expanded list file line that matches (iFile, SplitStep)
def find_ilist_and_total(master_list_dir, iFile, splitStep):
    """
    Return (ilist, SplitTotal) by scanning master lists under master_list_dir.
    Each master list line format is:
        ilist iFile SplitStep SplitTotal
    """
    if not os.path.isdir(master_list_dir):
        raise FileNotFoundError(f"Master list dir not found: {master_list_dir}")

    for fname in os.listdir(master_list_dir):
        if not fname.endswith('_list.list'):
            continue
        path = os.path.join(master_list_dir, fname)
        with open(path, 'r') as ml:
            for line in ml:
                parts = line.strip().split()
                if len(parts) < 4:
                    continue
                try:
                    file_ifile = int(parts[1])
                    file_step = int(parts[2])
                except ValueError:
                    continue
                if file_ifile == iFile and file_step == splitStep:
                    ilist = parts[0]
                    split_total = int(parts[3])
                    return ilist, split_total
    # Not found in expanded master lists: try to guess ilist via config/list filename
    # and conservatively return None so caller can decide.
    # attempt to find a per-root list file that ends with f'_{iFile}.list'
    if os.path.isdir(config_list_dir):
        for fname in os.listdir(config_list_dir):
            if fname.endswith(f'_{iFile}.list'):
                ilist_guess = os.path.join('config', 'list', DataSetName, fname)
                # We do not know split_total from this file alone - caller should handle
                return ilist_guess, None
    raise ValueError(f"Could not find ilist+SplitTotal for iFile={iFile}, splitStep={splitStep} in {master_list_dir}")

def makeSubmitScript(tuple_pairs, submitName, resubmit, maxResub, DataSetName):
    """
    Create a tuple submit file for resubmitting failed (iFile, SplitStep) jobs.
    Writes tuple_filelist with columns:
       ilist iFile SplitStep SplitTotal
    Then copies submitName_single.submit -> submitName_tuple.submit and replaces
    placeholders so the new submit reads the tuple file and submits the exact chunks.
    """
    tuple_pairs = sorted(tuple_pairs)
    tuple_filelist = f"{submitName}_tuple.txt"
    resubmitFiles = 0
    if tuple_pairs == []:
        return resubmitFiles  # nothing to do

    # derive workingDir from submitName (submitName == workingDir+"/src/"+DataSetName)
    workingDir = os.path.dirname(os.path.dirname(submitName))
    # location where the dataset master lists should live
    master_list_dir = os.path.join(workingDir, "list", DataSetName)
    # fallback config list location (per-root-list files)
    config_list_dir = os.path.join(workingDir, "config", "list", DataSetName)

    # Build tuple file containing ilist, iFile, SplitStep, SplitTotal
    with open(tuple_filelist, 'w') as tuple_file:
        for tuple_pair in tuple_pairs:
            resubmitFiles += 1
            if DO_EVENTCOUNT:
                tuple_file.write(f"{tuple_pair}\n")
                continue

            # tuple_pair assumed to be (iFile, splitStep)
            try:
                iFile, splitStep = tuple_pair
            except Exception:
                raise ValueError("tuple_pairs must contain (iFile, splitStep) pairs in non-eventcount mode")

            ilist, split_total = find_ilist_and_total(master_list_dir, iFile, splitStep)
            if split_total is None:
                # fallback for couldn't find SplitTotal, try to get it from the master lists by matching ilist
                # scan master lists for ilist match and take split_total from matching lines
                found = False
                for fname in os.listdir(master_list_dir):
                    if not fname.endswith('_list.list'):
                        continue
                    path = os.path.join(master_list_dir, fname)
                    with open(path, 'r') as ml:
                        for line in ml:
                            parts = line.strip().split()
                            if len(parts) >= 4 and parts[0].endswith(ilist.split('/')[-1]):
                                try:
                                    split_total = int(parts[3])
                                    found = True
                                    break
                                except Exception:
                                    continue
                    if found:
                        break
                if split_total is None:
                    raise ValueError(f"Could not determine SplitTotal for ilist {ilist}, iFile {iFile}, splitStep {splitStep}")

            # write the 4-column line: ilist iFile SplitStep SplitTotal
            tuple_file.write(f"{ilist} {iFile} {splitStep} {split_total}\n")

    # Prepare the tuple submit based on the single template
    newFileName = f"{submitName}_tuple.submit"
    os.system(f"cp {submitName}_single.submit {newFileName}")

    # Read the single submit and do safe, robust replacements
    with open(newFileName, 'r') as fh:
        file_content = fh.read()

    if DO_EVENTCOUNT:
        file_content = file_content.replace("0.list", "$(list).list")
        file_content = file_content.replace("0.root", "$(list).root")
        file_content = file_content.replace("0.out", "$(list).out")
        file_content = file_content.replace("0.err", "$(list).err")
        file_content = file_content.replace("0.log", "$(list).log")
        file_content = file_content.replace("queue", f"queue list from {tuple_filelist}")
    else:
        # replace the hard-coded single ifile path (the -ilist=...) with $(ilist)
        m = re.search(r'-ilist=([^\s]+)', file_content)
        if m:
            single_ifile = m.group(1)
            file_content = file_content.replace(single_ifile, '$(ilist)')
        else:
            # fallback: nothing found; assume single template used './config/list/..._0.list' literal 'X_0.list'
            file_content = file_content.replace("X_0.list", "$(ilist)")

        # 2) replace the literal 0_0 pattern in filenames with macro form $(iFile)_$(SplitStep)
        #    and ensure any other occurrences (output/error/log/transfer remap) follow suite
        file_content = file_content.replace("0_0", "$(iFile)_$(SplitStep)")
        # also handle older placeholder 0 (if present)
        file_content = file_content.replace("_0.list", "_$(iFile).list")

        # 3) replace the -split=1,<N> single-job spec with the macroized split spec
        #    find occurrences like '-split=1,123' and replace with '-split=$$([$(SplitStep)+1]),$(SplitTotal)'
        file_content = re.sub(r'-split=1,\s*\d+', '-split=$$([$(SplitStep)+1]),$(SplitTotal)', file_content)

        # 4) replace the queue line to include all four variables reading from tuple_filelist
        file_content = re.sub(r'\bqueue\b.*', f'queue ilist, iFile, SplitStep, SplitTotal from {tuple_filelist}', file_content, count=1)

    file_content = file_content.replace('RequestCpus = 1', 'RequestCpus = 2')
    file_content = file_content.replace('priority = 10', 'priority = 20')

    # write the new tuple submit
    with open(newFileName, 'w') as fh:
        fh.write(file_content)

    # Print and optionally submit
    print('total resubmit:', resubmitFiles, flush=True)
    if resubmit:
        if resubmitFiles > maxResub:
            print(f"You are about to make {resubmitFiles} and resubmit {resubmitFiles} jobs for dataset: {DataSetName}!", flush=True)
            print(f"You should double check there are no issues with your condor submissions.", flush=True)
            print(f"If you are confident you want to resubmit, then you should rerun this script with '-l {resubmitFiles}'.", flush=True)
            print(f"Or run condor_submit {newFileName}", flush=True)
        else:
            condor_monitor = CondorJobCountMonitor(threshold=0.99*get_auto_THRESHOLD()+resubmitFiles, verbose=False)
            condor_monitor.wait_until_jobs_below()
            os.system(f"condor_submit {newFileName}")

    return resubmitFiles

def testRootFile(root_file):
    status = True
    try:
        # Can ignore any message like: ReadStreamerInfo, class:string, illegal uid=-2
        # https://root-forum.cern.ch/t/readstreamerinfo-illegal-uid-with-newer-root-versions/41073
        ROOT.gErrorIgnoreLevel = ROOT.kFatal+1
        root_file_test = ROOT.TFile.Open(root_file);
        if not root_file_test or root_file_test.IsZombie() or root_file_test.GetListOfKeys().GetSize() == 0:
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

def ReadFilterList(filterlist_filename):
    try:
        with open(filterlist_filename, 'r') as filterlist:
            filterlist = [line.strip() for line in filterlist if line.strip()]
            if filterlist:
                return filterlist
            else:
                return []
    except FileNotFoundError:
        return None

def UpdateFilterList(DataSetName, filterlist_filename, add):
    # add == True means add to list, False means remove if in list
    try:
        with open(filterlist_filename,'r') as filterlist:
            dataset_list = {line.strip() for line in filterlist}
    except FileNotFoundError:
        dataset_list = set()  # If file doesn't exist, start with an empty set
    if add:
        dataset_list.add(DataSetName)  # Add dataset if not already present
    else:
        dataset_list.discard(DataSetName)  # Remove dataset if present
    # Write updated list back to the file
    with open(filterlist_filename, 'w') as filterlist:
        filterlist.write("\n".join(sorted(dataset_list)) + "\n")  # Sort for consistency

# Check condor jobs
def checkJobs(workingDir, outputDir, skipEC, skipDAS, skipMissing, skipSmall,
              skipErr, skipOut, skipZombie, resubmit, maxResub,
              filter_list, skipDASDataset, skipHists, doRucio):
    print("Running over the directory '{0}'.".format(workingDir), flush=True)
    print("------------------------------------------------------------", flush=True)
    grep_ignore = "-e \"Warning\" -e \"WARNING\" -e \"TTree::SetBranchStatus\" -e \"libXrdSecztn.so\" -e \"Phi_mpi_pi\" -e \"tar: stdout: write error\" -e \"INFO\""
    grep_ignore += " -e \"SetBranchAddress\""

    srcDir = os.listdir(os.path.join(workingDir, "src"))
    total_resubmit = 0

    # helper to parse tuples from a file path (base filename without extension)
    def path_to_tuple(path):
        """Return Tuple or int depending on DO_EVENTCOUNT. Safe parsing of basename."""
        base = os.path.basename(path)
        # remove common extensions if present
        for ext in (".err", ".out", ".root", ".list"):
            if base.endswith(ext):
                base = base[: -len(ext)]
        parts = base.split("_")
        try:
            if DO_EVENTCOUNT:
                return int(parts[-1])
            else:
                return (int(parts[-2]), int(parts[-1]))
        except Exception:
            # If parsing fails, return None so caller can ignore it.
            return None

    for file in srcDir:
        if "X.submit" not in file:
            continue
        # filter_list handling
        if filter_list:
            do_name = any(name in file for name in filter_list)
            if not do_name:
                continue
        DataSetName = file.split(".")[0]
        resubmit_set = set()
        event_count = EventCount()

        # --- DAS check ---
        if not skipDAS:
            nDAS_fail = 0
            NDAS_true = 0
            NDAS = 0
            listDir = os.path.join(workingDir, "config", "list", DataSetName)
            dataset = []
            try:
                for lFile in os.listdir(listDir):
                    if lFile.endswith('.list') and not lFile.endswith('_list.list'):
                        with open(os.path.join(listDir, lFile), 'r') as input_root_files:
                            input_root_filename = input_root_files.readline().strip()
                            dataset = event_count.GetDatasetFromFile(input_root_filename)
                            if dataset == []:
                                NDAS_true += event_count.EventsInDAS(input_root_filename, True)
                            else:
                                NDAS_true = event_count.EventsInDAS(dataset, False)
                                break
            except Exception:
                dataset = []

            # sum counts from outputs
            outfiles_dir = os.path.join(outputDir, DataSetName)
            try:
                outfiles = os.listdir(outfiles_dir)
                for outfile in outfiles:
                    NDAS += event_count.GetDASCount(os.path.join(outfiles_dir, outfile))
            except Exception:
                NDAS = 0

            if NDAS != NDAS_true:
                comp_percent = 0
                if NDAS_true > 0:
                    comp_percent = 100. * math.floor(NDAS / NDAS_true * 1000) / 1000
                print(f'{DataSetName} failed the DAS check! ({comp_percent}%)\nMissing {NDAS_true-NDAS} events', flush=True)
                if (not skipDASDataset) and dataset:
                    files_datasets = set()
                    try:
                        bash = f"find {os.path.join(outputDir, DataSetName)} -type f"
                        files = subprocess.check_output(['bash', '-c', bash], text=True).splitlines()
                        for f in files:
                            if not f:
                                continue
                            file_datasets = event_count.getFullDASDatasetNames(f)
                            for ds in file_datasets:
                                if "jsonparser" in ds:
                                    tup = path_to_tuple(f)
                                    if tup is None:
                                        continue
                                    if tup not in resubmit_set:
                                        resubmit_set.add(tup)
                                        nDAS_fail += 1
                                files_datasets.add(ds)
                    except Exception:
                        files_datasets = set()

                    dataset_set = set(dataset)
                    if files_datasets != dataset_set:
                        print(f'{DataSetName} dataset mismatch!', flush=True)
                        print(f'From files:   {sorted(files_datasets)}', flush=True)
                        print(f'From DAS: {sorted(dataset_set)}', flush=True)

                UpdateFilterList(DataSetName, f"{workingDir}/CheckFiles_FilterList.txt", True)
            else:
                print(f'{DataSetName} passed the DAS check!', flush=True)
                UpdateFilterList(DataSetName, f"{workingDir}/CheckFiles_FilterList.txt", False)
                # --- Histograms folder check ---
                # Need to do this regardless of DAS because it looks at different content in the root file
                num_hist = 0
                if not skipHists:
                    try:
                        bash = f"find {os.path.join(outputDir, DataSetName)} -type f"
                        histFiles = subprocess.check_output(['bash', '-c', bash], text=True).splitlines()
                    except subprocess.CalledProcessError:
                        histFiles = []
                    except Exception:
                        histFiles = []
                    for histFile in histFiles:
                        if not histFile:
                            continue
                        try:
                            if event_count.check_histograms_in_file(histFile) is True:
                                continue
                        except Exception:
                            # If histogram fails, treat as zombie and attempt to remove
                            pass
                        # attempt to remove the bad file
                        try:
                            os.remove(histFile)
                        except Exception:
                            pass
                        tup = path_to_tuple(histFile)
                        if tup is None:
                            continue
                        if tup not in resubmit_set:
                            resubmit_set.add(tup)
                            num_hist += 1
                if num_hist > 0:
                    print(f"Got {num_hist} root files fail histogram check for dataset {DataSetName}", flush=True)
                else:
                    continue # passed DAS and hist checks -> skip other checks
            print(f"Got {nDAS_fail} DAS fail files for dataset {DataSetName}", flush=True)

        # --- missing files check ---
        if not skipMissing:
            listDir = os.path.abspath(os.path.join(workingDir, "list", DataSetName))
            listFiles = []
            try:
                listFiles = os.listdir(listDir)
            except Exception:
                listFiles = []
            added_missing = 0
            if DO_EVENTCOUNT:
                numList = len([lf for lf in listFiles if os.path.isfile(os.path.join(listDir, lf))]) - 1
                missing = getMissingFilesEC(outputDir+'/'+DataSetName, numList)
                for t in missing:
                    if t not in resubmit_set:
                        resubmit_set.add(t)
                        added_missing += 1
            else:
                try:
                    missing = getMissingFiles(os.path.join(outputDir, DataSetName), os.path.abspath(os.path.join(workingDir, "list", DataSetName, DataSetName + "_list.list")))
                except Exception:
                    missing = []
                for t in missing:
                    if t not in resubmit_set:
                        resubmit_set.add(t)
                        added_missing += 1
            print(f"Got {added_missing} missing files for dataset {DataSetName}", flush=True)

        # --- small files check ---
        if not skipSmall:
            bash = f"find {os.path.join(outputDir, DataSetName)} -type f -size -1k"
            try:
                smallFiles_out = subprocess.check_output(['bash', '-c', bash], text=True).splitlines()
            except subprocess.CalledProcessError:
                smallFiles_out = []
            except Exception:
                smallFiles_out = []

            num_small = 0
            for smallFile in smallFiles_out:
                if not smallFile:
                    continue
                tup = path_to_tuple(smallFile)
                if tup is None:
                    continue
                if tup not in resubmit_set:
                    resubmit_set.add(tup)
                    num_small += 1
            print(f"Got {num_small} small files for dataset {DataSetName}", flush=True)

        # --- .err files check ---
        if not skipErr:
            marker = "WARNING: In non-interactive mode release checks e.g. deprecated releases, production architectures are disabled."
            num_error = 0

            for digit in range(10):
                pattern = os.path.join(workingDir, "err", DataSetName, f"*{digit}.err")
                files = glob.glob(pattern)  # <-- expand the glob in Python
            
                if not files:
                    continue  # skip if no matching files
            
                # Join files into a single space-separated string
                file_list = " ".join(files)
            
                awk_cmd = (
                    f"awk -v marker='{marker}' "
                    "'{ if (found) { print FILENAME \":\" $0 } } "
                    "$0 ~ marker { found=1 }' "
                    f"{file_list}"
                )
            
                bash = f"{awk_cmd} | grep -v {grep_ignore}"
            
                try:
                    output = subprocess.check_output(['bash', '-c', bash], text=True).splitlines()
                except subprocess.CalledProcessError:
                    output = []
                except Exception:
                    output = []
            
                for line in output:
                    if not line:
                        continue
                    filename_part = line.split(":", 1)[0]
                    tup = path_to_tuple(filename_part)
                    if tup is None:
                        continue
                    if tup not in resubmit_set:
                        resubmit_set.add(tup)
                        num_error += 1

            print(f"Got {num_error} error files for dataset {DataSetName}", flush=True)

        # --- .out files check ---
        if not skipOut:
            bash = f"grep -r \"Ntree 0\" {os.path.join(workingDir, 'out', DataSetName)}"
            try:
                outLines = subprocess.check_output(['bash', '-c', bash], text=True).splitlines()
            except subprocess.CalledProcessError:
                outLines = []
            except Exception:
                outLines = []
            num_out = 0
            for outLine in outLines:
                if not outLine:
                    continue
                filename_part = outLine.split(":", 1)[0]
                tup = path_to_tuple(filename_part)
                if tup is None:
                    continue
                if tup not in resubmit_set:
                    resubmit_set.add(tup)
                    num_out += 1
            print(f"Got {num_out} out files for dataset {DataSetName}", flush=True)

        # --- zombie root files check ---
        if not skipZombie:
            num_zomb = 0
            try:
                bash = f"find {os.path.join(outputDir, DataSetName)} -type f"
                zombFiles = subprocess.check_output(['bash', '-c', bash], text=True).splitlines()
            except subprocess.CalledProcessError:
                zombFiles = []
            except Exception:
                zombFiles = []
            for zombFile in zombFiles:
                if not zombFile:
                    continue
                try:
                    if testRootFile(zombFile) is True:
                        continue
                except Exception:
                    # If testRootFile fails, treat as zombie and attempt to remove
                    pass
                # attempt to remove the bad file
                try:
                    os.remove(zombFile)
                except Exception:
                    pass
                tup = path_to_tuple(zombFile)
                if tup is None:
                    continue
                if tup not in resubmit_set:
                    resubmit_set.add(tup)
                    num_zomb += 1
            print(f"Got {num_zomb} zombie root files for dataset {DataSetName}", flush=True)

        # --- EventCount check (EC) ---
        if not skipEC:
            try:
                failedEC = eventCountCheck(os.path.join(outputDir, DataSetName))
            except Exception:
                failedEC = []
            added_ec = 0
            for t in failedEC:
                if t not in resubmit_set:
                    resubmit_set.add(t)
                    added_ec += 1
            print(f"Got {len(failedEC)} files failed EventCount check for dataset {DataSetName}", flush=True)

        # create submit script(s) and count new jobs
        nJobs = makeSubmitScript(resubmit_set, os.path.join(workingDir, "src", DataSetName),
                                 resubmit, maxResub, DataSetName)

        # create Rucio request for bad files
        if doRucio:
            # Setup rucio environment once per dataset
            rucio_setup = "source /cvmfs/cms.cern.ch/rucio/setup-py3.sh"

            # Path to tuple file
            tuple_file = os.path.join(
                workingDir,
                "src",
                f"{DataSetName}_tuple.txt"
            )

            if not os.path.exists(tuple_file):
                print(f"[Rucio] Tuple file not found: {tuple_file}", flush=True)
            else:
                rucio_files = set()

                with open(tuple_file, "r") as tf:
                    for line in tf:
                        if not line.strip():
                            continue

                        parts = line.strip().split()
                        if len(parts) < 2:
                            continue

                        list_path = parts[0]

                        # Full path to the .list file
                        full_list_path = os.path.join(workingDir, list_path)

                        if not os.path.exists(full_list_path):
                            print(f"[Rucio] Missing list file: {full_list_path}", flush=True)
                            continue

                        with open(full_list_path, "r") as lf:
                            for rootline in lf:
                                rootline = rootline.strip()
                                if not rootline:
                                    continue

                                # Strip everything before /store/
                                store_idx = rootline.find("/store/")
                                if store_idx == -1:
                                    continue

                                rucio_path = "cms:" + rootline[store_idx:]
                                if has_rule(rucio_path, "T2_US_Nebraska"):
                                    continue
                                rucio_files.add(rucio_path)
                if len(rucio_files) > 0:
                    print(f"Creating {len(rucio_files)} rucio rules")
                    for rfile in sorted(rucio_files):
                        cmd = (
                            f'{rucio_setup} && '
                            f'rucio rule add --copies 1 --lifetime 175000 --rses T2_US_Nebraska -d {rfile}'
                        )
                        try:
                            subprocess.check_call(['bash', '-c', cmd])
                        except subprocess.CalledProcessError as e:
                            print(f"[Rucio] Failed rule for {rfile}", flush=True)

        total_resubmit += nJobs

    return total_resubmit

def main():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--eventCount",           action='store_true', help="for checking event count jobs")
    parser.add_argument("--directory",      "-d", default=None, help="working directory for condor submission (ex: Summer23BPix_130X)")
    parser.add_argument("--output",         "-o", default="/ospool/cms-user/"+USER+"/NTUPLES/Processing/Summer23BPix_130X/", help="output area for root files")
    parser.add_argument("--checkAll",       "-a", action='store_true', help="explicitly check all datasets in working directory (overwrite past checks)")
    parser.add_argument("--skipResubmit",   "-r", action='store_false', help="do not resubmit files")
    parser.add_argument("--skipEC",         "-c", action='store_true', help="skip checking event count tree against internal DAS")
    parser.add_argument("--skipDAS",        "-w", action='store_true', help="skip checking compared to central DAS")
    parser.add_argument("--skipMissing",    "-m", action='store_true', help="skip checking missing files")
    parser.add_argument("--skipSmall",      "-s", action='store_true', help="skip checking small files")
    parser.add_argument("--skipErr",        "-e", action='store_true', help="skip checking err files (recommended: slow)")
    parser.add_argument("--skipOut",        "-u", action='store_true', help="skip checking out files")
    parser.add_argument("--skipZombie",     "-z", action='store_true', help="skip checking zombie root files")
    parser.add_argument("--skipDASDataset", "-k", action='store_true', help="skip checking DAS dataset matches")
    parser.add_argument("--skipHists",      "-f", action='store_true', help="skip checking Histograms folder")
    parser.add_argument("--rucio",          "-n", action='store_true', help="create rucio request to move bad files to different site")
    parser.add_argument("--maxResub",       "-l", default=5000, help="max number of jobs to resubmit")
    parser.add_argument("--threshold",      "-t", default=0.99*get_auto_THRESHOLD(), help="min number of jobs running before starting checker")
    parser.add_argument("--sleep",          "-p", default=1, help="time to sleep before starting checker")

    global DO_EVENTCOUNT
    options        = parser.parse_args()
    DO_EVENTCOUNT  = options.eventCount
    directory      = options.directory	
    output         = options.output
    checkAll       = options.checkAll
    resubmit       = options.skipResubmit
    skipEC         = options.skipEC
    skipDAS        = options.skipDAS
    skipMissing    = options.skipMissing
    skipSmall      = options.skipSmall
    skipErr        = options.skipErr
    skipOut        = options.skipOut
    skipZombie     = options.skipZombie
    skipDASDataset = options.skipDASDataset
    skipHists      = options.skipHists
    doRucio        = options.rucio
    maxResub       = int(options.maxResub)
    threshold      = int(options.threshold)
    sleep_time     = int(options.sleep)

    # quick strictly das check from root file (useful for after final hadd test)
    if directory is None:
        file_list = []
        event_count = EventCount()

        if "root://cmseos.fnal.gov/" in output:
            eos_path = output.split("root://cmseos.fnal.gov/")[1]
            # Run xrdfs and capture stdout
            result = subprocess.run(
                ["xrdfs", "root://cmseos.fnal.gov", "ls", eos_path],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            file_list = [
                f"root://cmseos.fnal.gov/{line.strip()}"
                for line in result.stdout.splitlines()
                if line.strip().endswith(".root")
            ]
        else:
            # Local directory
            file_list = [
                os.path.join(output, f)
                for f in os.listdir(output)
                if f.endswith(".root") and os.path.isfile(os.path.join(output, f))
            ]

        for full_path in file_list:
            totalEvents = event_count.countTotalEvents(full_path)
            DAS_events = event_count.getEventsFromDASDatasetNames(full_path)
            if DAS_events == 0:
                print("Something wrong with DAS events in", full_path, flush=True)
            elif totalEvents != DAS_events:
                comp_percent = 100. * math.floor(totalEvents / DAS_events * 1000) / 1000
                print(f'{full_path} failed the DAS check! ({comp_percent}%) Use other options to investigate', flush=True)
            else:
                print(f'{full_path} passed the DAS check!', flush=True)

    else:
        # default checking that's robust
        filter_list = [ # list of name of samples to explicitly check (empty list will check all)
            #Example: "WZ_TuneCP5_13p6TeV_pythia8",
        ]
        fileFilterList = ReadFilterList(f"{directory}/CheckFiles_FilterList.txt")
        if fileFilterList is not None: filter_list.extend(fileFilterList)
        if fileFilterList == [] and not checkAll:
            print("All datasets passed DAS check!", flush=True)
            return
        if checkAll: filter_list.clear()

        time.sleep(sleep_time)
        condor_monitor = CondorJobCountMonitor(threshold=threshold,verbose=False)
        print(f"Waiting until minumum of {threshold} jobs in the queue", flush=True)
        condor_monitor.wait_until_jobs_below()
        print("Running checker...", flush=True)
        nJobs = checkJobs(directory,output,skipEC,skipDAS,skipMissing,skipSmall,skipErr,skipOut,skipZombie,resubmit,maxResub,filter_list,skipDASDataset,skipHists,doRucio)
        if resubmit and nJobs > 0:
            print(f"Resubmitted a total of {nJobs} jobs!", flush=True)
    print("Checker Complete!", flush=True)

if __name__ == "__main__":
    main()
