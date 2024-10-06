# tools.py

import os
import csv
import glob
import ROOT

# TODO

# DONE
# - Update get_eos_file_list() to use a pattern

# creates directory if it does not exist
def makeDir(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

# get numbers from a string using a delimiter separation
def getNumbers(string, delimiter):
    numbers = [int(s) for s in string.split(delimiter) if s.isdigit()]
    return numbers

# determine if number is in range
def numberInRange(number, range_min, range_max):
    if number >= range_min and number <= range_max:
        return True
    else:
        return False

# Note for csv library:
# The "newline" argument is available in python 3, but not in python 2.

# read csv file: takes a csv file as input and outputs data in a matrix
def readCSV(input_file):
    data = []
    with open(input_file, mode="r") as f:
        reader = csv.reader(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for row in reader:
            data.append(row)
    return data

# write csv file: takes data matrix as input and outputs a csv file 
def writeCSV(output_file, data):
    if '/' in output_file:
        csv_dir = output_file.split("/")[0:-1]
        csv_dir = '/'.join(csv_dir)
        os.system(f"mkdir -p {csv_dir}")
    with open(output_file, mode="w") as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for row in data:
            writer.writerow(row)

# get chain from list of ROOT files
def getChain(input_files, num_files):
    verbose = False
    
    # use num_files as max if it is not negative
    if num_files >= 0:
        input_files = input_files[0:num_files]
    n_input_files = len(input_files)
    
    # Create TChain
    chain = ROOT.TChain('pixelTree')
    
    # Add files to chain
    if verbose:
        print("Adding {0} file(s) to chain:".format(n_input_files))
    for f in input_files:
        if verbose:
            print(" - {0}".format(f))
        chain.Add(f)
    
    return chain

# get a list of subdirectories in a directory
def get_subdirectories(directory):
    subdirectories = [d for d in os.listdir(directory) if os.path.isdir(os.path.join(directory, d))]
    return subdirectories

# count files in a directory matching a pattern
def count_files_with_pattern(directory, pattern):
    files = glob.glob(os.path.join(directory, pattern))
    n_files = len(files)
    return n_files

# get list of local files
def get_file_list(dir_):
    file_list_tmp = [os.path.join(dir_, f) for f in os.listdir(dir_) if (os.path.isfile(os.path.join(dir_, f)) and ('.root' in f))]
    return file_list_tmp

# get list of local files using glob
def get_file_list_glob(directory, pattern="*.root"):
    d = directory
    if d[-1] != "/":
        d += "/"
    return glob.glob(d + pattern)

# get list of ROOT files on EOS
# - if pattern is set, require that pattern is in base name
def get_eos_file_list(path, pattern="", eosurl="root://cmseos.fnal.gov"):
    debug = False
    output = [] 
    with eosls(path, "", eosurl) as files:
        for f in files:
            name        = f.strip()
            full_name   = "{0}/{1}".format(eosurl, name)
            base_name   = os.path.basename(name) 
            
            # print info for debugging
            if debug:
                print("In get_eos_file_list():")
                print(" - name: {0}".format(name))
                print(" - full_name: {0}".format(full_name))
                print(" - base_name: {0}".format(base_name))
            
            # require file to be a ROOT file 
            if name.endswith(".root"):
                # if pattern is set, require that pattern is in base name
                if pattern:
                    if pattern in base_name:
                        output.append(full_name)
                else:
                    output.append(full_name)
    return output

# eosls command using xrdfs
def eosls(path, option="", eosurl="root://cmseos.fnal.gov"):
    debug = True
    command = "xrdfs %s ls %s %s" % (eosurl, option, path)
    if debug:
        print(command)
    return os.popen(command)

# eosrm command using xrdfs
def eosrm(path, option="", eosurl="root://cmseos.fnal.gov"):
    debug = True
    command = "xrdfs %s rm %s %s" % (eosurl, option, path)
    if debug:
        print(command)
    return os.popen(command)
