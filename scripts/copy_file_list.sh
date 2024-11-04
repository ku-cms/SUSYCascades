#!/usr/bin/bash
# copy_file_list.sh

# Copies all files in a file list to an output directory using xrdcp.

# Requires CERN grid certificate and active voms proxy:
# voms-proxy-init --valid 192:00 -voms cms
# voms-proxy-info

# To run interactively:
# ./copy_file_list.sh file_list.txt output_dir

# To run using nohup (useful when it takes a long time):
# nohup ./copy_file_list.sh file_list.txt output_dir > copy_list_001.log 2>&1 & 

# User inputs:
# - file list
# - output directory

# Output:
# - copies files to output directory

# User inputs:
FILE_LIST=$1
OUTPUT_DIR=$2

# Main FNAL redirector (searches all sites):
REDIRECTOR=root://cmsxrootd.fnal.gov/

START_TIME=$EPOCHREALTIME

echo "FILE_LIST: ${FILE_LIST}"
echo "OUTPUT_DIR: ${OUTPUT_DIR}"
echo "REDIRECTOR: ${REDIRECTOR}"

# Check if FILE_LIST is empty.
if [[ -z "$FILE_LIST" ]]; then
    echo "ERROR: FILE_LIST is empty."
    echo "Please provide a file list as the first argument."
    exit 1
fi

# Check if OUTPUT_DIR is empty.
if [[ -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: OUTPUT_DIR is empty."
    echo "Please provide an ouput directory as the second argument."
    exit 1
fi

# Check if file list does not exist.
if [[ ! -f "$FILE_LIST" ]]; then
    echo "ERROR: The file '${FILE_LIST}' does not exist!"
    exit 1
fi

# Create output directory.
mkdir -p ${OUTPUT_DIR}

echo "Copying files using xrdcp..."

while read line; do
    echo " - ${line}"
    xrdcp -f ${REDIRECTOR}${line} ${OUTPUT_DIR}
done < ${FILE_LIST}

# Print number of files copied
NUM_FILES=$(ls ${OUTPUT_DIR} | wc -l)
echo "Number of files in output directory: ${NUM_FILES}"

# Print storage size of output directory
OUTPUT_SIZE=$(du -sh ${OUTPUT_DIR})
echo "Storage size of output directory: ${OUTPUT_SIZE}"

# Calculate run time
END_TIME=$EPOCHREALTIME
RUN_TIME_SEC=$(bc -l <<< "$END_TIME - $START_TIME")
RUN_TIME_MIN=$(bc -l <<< "$RUN_TIME_SEC / 60")
RUN_TIME_HR=$(bc -l <<< "$RUN_TIME_MIN / 60")

# Print run time
printf "Run time: %0.2f seconds = %0.2f minutes = %0.2f hours\n" ${RUN_TIME_SEC} ${RUN_TIME_MIN} ${RUN_TIME_HR}

echo "Done!"

