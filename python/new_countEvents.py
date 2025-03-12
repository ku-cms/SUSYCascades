# countEvents.py

import os
import glob
import ROOT
import argparse
from colorama import Fore, Back, Style
import tools
import subprocess
import json

# Count the number of "total" and "saved" events for all ROOT files in a directory.
#
# Example for DAS check on connect:
#
# python3 python/new_countEvents.py -d /ospool/cms-user/zflowers/NTUPLES/Processing/Summer23BPix_130X/ -w
#
# Check event count file:
#
# python3 python/new_countEvents.py --eventCount --filetag Summer23_130X
# 
# This script is intended for checking counts after the fact and only does "internal" DAS counts
# For running on newly processed jobs, recommend using python/CheckFiles.py
# CheckFiles is more robust and also has option for "offline" DAS checking
# Functions called here only use "online" DAS check which looks at what's saved in the output root file

# Make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# Make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)
# Tell ROOT not to be in charge of memory, fix issue of histograms being deleted when ROOT file is closed:
ROOT.TH1.AddDirectory(False)

class EventCount:
    def __init__(self, event_count_tree="EventCount", analysis_tree="KUAnalysis", base_event_count=""):
        self.event_count_tree   = event_count_tree
        self.analysis_tree      = analysis_tree

        self.base_event_count   = base_event_count
        if base_event_count != "":
            root_file = "root/EventCount/EventCount_NANO_"+self.base_event_count
            root_file += ".root"
            self.base_analysis_tree_map = self.LoadEventCountMap(root_file)

    def LoadEventCountMap(self, root_file):
        try:
            # Can ignore any message like: ReadStreamerInfo, class:string, illegal uid=-2
            # https://root-forum.cern.ch/t/readstreamerinfo-illegal-uid-with-newer-root-versions/41073
            root_file_test = ROOT.TFile.Open(root_file);
            if root_file_test.IsOpen():
                root_file_test.Close()
            if not root_file_test or root_file_test.IsZombie():
                return None
        except:
            return None
        analysis_tree_map = {}
        tree = self.GetEventCountTree()
        chain = ROOT.TChain(tree)
        chain.Add(root_file)
        n_entries = chain.GetEntries()
        for i in range(n_entries):
            chain.GetEntry(i)
            keyList = []
            keyList.append(chain.dataset)
            keyList.append(str(chain.MP))
            keyList.append(str(chain.MC))
            key = ''.join(k+"_" for k in keyList)
            key = key[:-1] 
            if key in analysis_tree_map.keys():
                analysis_tree_map[key] += chain.Nevent
            else:
                analysis_tree_map[key] = chain.Nevent
        return analysis_tree_map

    def GetEventCountTree(self):
        return self.event_count_tree
    
    def SetEventCountTree(self, event_count_tree):
        self.event_count_tree = event_count_tree
    
    def GetAnalysisTree(self):
        return self.analysis_tree
    
    def SetAnalysisTree(self, analysis_tree):
        self.analysis_tree = analysis_tree
    
    # Gets "on-the-fly" DAS count from root file
    def GetDASCount(self, root_file):
        DAS_count = 0
        tree = self.GetEventCountTree()
        chain = ROOT.TChain(tree)
        chain.Add(root_file)
        n_entries = chain.GetEntries()
        for i in range(n_entries):
            chain.GetEntry(i)
            n_DAS_Count = chain.NDAS
            DAS_count += n_DAS_Count
        return int(DAS_count)

    def checkInternalDASCount(self, file):
        NDAS = self.GetDASCount(file)   
        Nevent = self.countTotalEvents(file)
        return NDAS == Nevent

    # Gets events directly from DAS (slow)
    def EventsInDASFile(self, u_file):
        """Get the number of events in DAS for the given file."""
        events = 0
        pos = u_file.find("/store/")
        if pos != -1:
            filename = u_file[pos:]  # Remove everything before "/store/"
        else:
            filename = u_file
        try:
            das_output = subprocess.check_output(
                ["dasgoclient", "-query", f'file={filename}', "-json"],
                text=True
            )
            das_data = json.loads(das_output)
            # Extract number of events from JSON
            if "nevents" in das_output:
                events = int(das_output.split('"nevents":', 1)[1].split(',', 1)[0])
        except (subprocess.CalledProcessError, ValueError, json.JSONDecodeError) as e:
            print(f"Error querying DAS: {e}")
        return events

    # Gets events directly from DAS
    def EventsInDASDataset(self, u_dataset):
        """Get the number of events in DAS for the given dataset."""
        events = 0
        try:
            query = f'dataset={u_dataset}'
            das_output = subprocess.check_output(
                ["dasgoclient", "-query", query, "-json"],
                text=True
            )
            das_data = json.loads(das_output)
            # Extract number of events from JSON
            if "nevents" in das_output:
                events = int(das_output.split('"nevents":', 1)[1].split(',', 1)[0])
        except (subprocess.CalledProcessError, ValueError, json.JSONDecodeError) as e:
            print(f"Error querying DAS: {e}")
        return events

    def EventsInDASDatasets(self, u_datasets):
        """Get the number of events in DAS for the given datasets."""
        events = 0
        datasets = u_datasets
        if type(datasets) is not list:
            datasets = [datasets]
        for dataset in datasets:
            events += self.EventsInDASDataset(dataset)
        return events

    def GetDatasetFromFile(self, u_file):
        """Get dataset name from file"""
        pos = u_file.find("/store/")
        if pos != -1:
            filename = u_file[pos:]  # Remove everything before "/store/"
        else:
            filename = u_file
        try:
            das_output = subprocess.check_output(
                ["dasgoclient", "-query", f'dataset file={filename}'],
                text=True,
                stderr=subprocess.STDOUT
            ).strip()
            # need to search more generically to check for exts
            das_output = das_output.split('/')
            dataset_name = das_output[1]
            campaign_tags = das_output[2]
            campaign_tags = campaign_tags.split('-')[0]
            aod_version = das_output[3]
            query = f'dataset=/{dataset_name}/{campaign_tags}*/{aod_version}'
            das_output = subprocess.check_output(
                ["dasgoclient", "-query", query],
                text=True,
                stderr=subprocess.STDOUT
            ).strip()
            das_output = das_output.split('\n')
            for dataset in das_output:
                if 'JME' in dataset or 'PUFor' in dataset:
                    das_output.remove(dataset)
            return das_output
        except subprocess.CalledProcessError as e:
            print(f"Error querying DAS: {e.output.strip()}")
            return []

    # Assume user is passing file, user can override to check for entire dataset
    def EventsInDAS(self, u_input, file=True):
        """Calls needed helper function based on file bool"""
        if file:
            return self.EventsInDASFile(u_input)
        else:
            return self.EventsInDASDatasets(u_input)    

    # count total events in a ROOT file
    # iterate over entries in the event count tree
    def countTotalEvents(self, root_file):
        result = 0
        tree = self.GetEventCountTree()
        chain = ROOT.TChain(tree)
        chain.Add(root_file)
        n_entries = chain.GetEntries()
        for i in range(n_entries):
            chain.GetEntry(i)
            n_events = chain.Nevent
            result += n_events
        return int(result)

    # count saved events in a ROOT file
    # use the number of entries in the analysis tree
    def countSavedEvents(self, root_file):
        tree = self.GetAnalysisTree()
        chain = ROOT.TChain(tree)
        chain.Add(root_file)
        n_events = chain.GetEntries()
        return int(n_events)

    def checkEventCountFile(self, filetag):
        root_file = f"root/EventCount/EventCount_NANO_{filetag}.root"
        if not os.path.isfile(root_file):
            print(f"file: {root_file} does not exist!")
            return
        tree = self.GetEventCountTree()
        chain = ROOT.TChain(tree)
        chain.Add(root_file)
        n_entries = chain.GetEntries()
        dataset_counts = {}
        DAS_counts = {}
        for i in range(n_entries):
            chain.GetEntry(i)
            dataset = str(chain.dataset)
            Nevent = int(chain.Nevent)
            try:
                NDAS = int(chain.NDAS)
            except:
                print(f'input root file does not have DAS values saved!') 
                return
            if dataset in dataset_counts:
                dataset_counts[dataset] += Nevent
                DAS_counts[dataset] += NDAS
            else:
                dataset_counts[dataset] = Nevent
                DAS_counts[dataset] = NDAS
        for dataset in DAS_counts:
            if DAS_counts[dataset] != dataset_counts[dataset]:
                if DAS_counts[dataset] > 0:
                    print(Fore.RED + f"dataset: {dataset} has failed the check!" + Fore.RESET)
                    perc = round(100.*dataset_counts[dataset]/DAS_counts[dataset],2)
                    print(f"dataset: {dataset} is at {perc}%")
                else:
                    print(f"Got {DAS_counts[dataset]} events from the DAS check for {dataset}")
                    print(Fore.RED + f"dataset: {dataset} has failed the check!" + Fore.RESET)
            else:
                print(Fore.GREEN + f"dataset: {dataset} passes the DAS check" + Fore.RESET)

    # process directory containing ROOT files
    def processDir(self, directory, pattern, csv, sms, eos, verbose, das):
        if verbose:
            print(Fore.GREEN + "Counting events." + Fore.RESET)
            print("----------------------------")
            print("directory: {0}".format(directory))
            print("pattern: {0}".format(pattern))
            print("csv: {0}".format(csv))
            print("sms: {0}".format(sms))
            print("eos: {0}".format(eos))
            print("verbose: {0}".format(verbose))
            print("das: {0}".format(das))
            print("----------------------------")
        
        root_files = []
        base_file_names = []
        output_data = []
        n_events_map = {}
        sum_total_events = 0
        sum_saved_events = 0
        sum_das_events   = 0
        das_events       = 0
        
        # get ROOT files
        # - if pattern is set, then require file name to contain pattern
        if eos:
            root_files = tools.get_eos_file_list(directory, pattern)
        else:
            if pattern:
                root_files = glob.glob("{0}/*{1}*.root".format(directory, pattern))
            else:
                root_files = glob.glob("{0}/*.root".format(directory))

        # sort ROOT files alphabetically    
        root_files.sort()
        n_root_files = len(root_files)

        if verbose:
            #print("ROOT files: {0}".format(root_files))
            print(Fore.GREEN + "Found {0} ROOT files with these base names:".format(n_root_files) + Fore.RESET)

        # count events
        ntuple_analysis_tree_map = {}
        first_rt_file = root_files[0]
        for root_file in root_files:
            file_analysis_tree_map = self.LoadEventCountMap(root_file)
            if file_analysis_tree_map is None:
                continue
            for key,value in file_analysis_tree_map.items():
                if key in ntuple_analysis_tree_map:
                    ntuple_analysis_tree_map[key] += value
                else:
                    ntuple_analysis_tree_map[key] = value
            base_name = os.path.basename(root_file)
            base_name = base_name.replace(".root", "")
            base_file_names.append(base_name)
            n_saved_events = 0
            n_total_events = 0
            if sms:
                for key in file_analysis_tree_map:
                    tree = key.split("_")
                    tree = tree[-2:] # just mass points
                    tree = "".join(t+"_" for t in tree)
                    tree = "SMS_"+tree
                    tree = tree[:-1]
                    self.SetAnalysisTree(tree)
                    n_file_saved_events = self.countSavedEvents(root_file)        
                    n_saved_events += n_file_saved_events
            else:
                n_saved_events = self.countSavedEvents(root_file)
            # Get DAS Count
            if(das):
                das_events += self.GetDASCount(root_file)
            n_total_events = self.countTotalEvents(root_file)
            n_events_map[base_name] = {}
            n_events_map[base_name]["n_total_events"] = n_total_events
            n_events_map[base_name]["n_saved_events"] = n_saved_events
            if sms:
                n_events_map[base_name]["analysis_tree"] = self.analysis_tree
            sum_total_events += n_total_events
            sum_saved_events += n_saved_events
            if verbose:
                print(" - {0}".format(base_name))

        for key in ntuple_analysis_tree_map:
            if self.base_event_count == "":
                root_file = first_rt_file
                try:
                    # Can ignore any message like: ReadStreamerInfo, class:string, illegal uid=-2
                    # https://root-forum.cern.ch/t/readstreamerinfo-illegal-uid-with-newer-root-versions/41073
                    root_file_test = ROOT.TFile.Open(root_file);
                    if root_file_test.IsOpen():
                        root_file_test.Close()
                    if not root_file_test or root_file_test.IsZombie():
                        print("can't open first root_file {root_file} to get filetag!")
                except:
                    print("can't open first root_file {root_file} to get filetag!")
                tree = "EventCount"
                chain = ROOT.TChain(tree)
                chain.Add(root_file)
                chain.GetEntry(0)
                self.base_event_count = str(chain.filetag)
                ec_root_file = "root/EventCount/EventCount_NANO_"+self.base_event_count
                if sms:
                    ec_root_file += "_SMS"
                ec_root_file += ".root"
                self.base_analysis_tree_map = self.LoadEventCountMap(ec_root_file)
            if key in self.base_analysis_tree_map:
                row = [self.base_event_count,key]
                base_events = self.base_analysis_tree_map[key]
                ntuple_saved_events = ntuple_analysis_tree_map[key]
                row.append(round(ntuple_saved_events/base_events,4))
                if base_events != ntuple_saved_events and verbose:
                    print(f"{key} missing {ntuple_saved_events/base_events} events")
                output_data.append(row)
        
        # print results
        for base_name in base_file_names:
            analysis_tree   = ""
            n_total_events = n_events_map[base_name]["n_total_events"]
            n_saved_events = n_events_map[base_name]["n_saved_events"]
            if verbose:
                if sms:
                    analysis_tree = n_events_map[base_name]["analysis_tree"]
                if sms:
                    print("{0}: {1}, DAS: {2}, {3}: {4}".format(base_name, n_total_events, das_events, analysis_tree, n_saved_events))
                else:
                    print("{0}: {1}, DAS: {2}, Saved: {3}".format(base_name, n_total_events, das_events, n_saved_events))
        if(das):
            base_name = base_file_names[0]
            das_percent = 100.*sum_total_events/das_events
            dataset = base_name
            if(das_events == sum_total_events):
                print(f"{dataset} passed the DAS check")
            else:
                print(f"{dataset} failed the DAS check and only processed {das_percent}% of the total dataset!")
        print("Sum of total events from all samples: {0}".format(sum_total_events))
        print("Sum of saved events from all samples: {0}".format(sum_saved_events))
        
        # if csv file name is set, then save data to csv file
        if csv:
            tools.writeCSV(csv, output_data)

def run():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--directory",  "-d", default="",                               help="directory containing ROOT files")
    parser.add_argument("--pattern",    "-p", default="",                               help="pattern for root file names (optional)")
    parser.add_argument("--csv",        "-c", default="",                               help="output csv file name (optional)")
    parser.add_argument("--sms",        "-s", default = False,  action = "store_true",  help="run over signal sample (optional)")
    parser.add_argument("--eos",        "-e", default = False,  action = "store_true",  help="run over ROOT files on EOS")
    parser.add_argument("--verbose",    "-v", default = False,  action = "store_true",  help="verbose flag to print more things")
    parser.add_argument("--das",        "-w", default = False,  action = "store_true",  help="get DAS count")
    parser.add_argument("--eventCount", "-t", default = False,  action = "store_true",  help="check event count file")
    parser.add_argument("--filetag",    "-f", default = "",                             help="filetag for event count checking (required for --eventCount check)")

    options     = parser.parse_args()
    directory   = options.directory
    pattern     = options.pattern
    csv         = options.csv
    sms         = options.sms
    eos         = options.eos
    verbose     = options.verbose
    das         = options.das
    eventCount  = options.eventCount
    filetag     = options.filetag

    if eventCount:
        event_count = EventCount()
        if directory: # option to pass directory to check intermediate event count root files
            for dirpath, _, files in os.walk(directory):
                for file in files:
                    if file.endswith('.root'):
                        if event_count.checkInternalDASCount(os.path.join(dirpath, file)): print(os.path.join(dirpath, file),"fails the DAS check!")
        # Check that the filetag is set.
        elif not filetag:
            print("You asked to check event count (--eventCount) but did not provide the required option filetag (--filetag FILETAG).")
            print("Please provide a filetag, for example: --filetag Summer23_130X")
            return
        else:
            event_count.checkEventCountFile(filetag)

    else:
        # Check that the directory is set.
        if not directory:
            print(Fore.RED + "ERROR: 'directory' is not set. Please provide a directory using the -d option." + Fore.RESET)
            return
        
        if sms:
            event_count = EventCount(event_count_tree="EventCount", analysis_tree="", base_event_count=filetag)
        else:
            event_count = EventCount(base_event_count=filetag)
        # Check if directory has root files. If no root files, walk down to root files and loop over subdirs
        # Need to have this just stop at top level directory with root files (no need to check staged hadded files (*_0.root, *_1.root, etc.)
        for dirpath, _, files in os.walk(directory):
            if any(f.endswith('.root') for f in files):
                event_count.processDir(dirpath, pattern, csv, sms, eos, verbose, das)
                break # Stop walking down once founded last stage of hadd files(?)

if __name__ == "__main__":
    run()

