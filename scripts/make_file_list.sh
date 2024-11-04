#!/usr/bin/bash
# make_file_list.sh

# Runs dasgoclient query for a given dataset to create a file list.

# Requires CERN grid certificate and active voms proxy:
# voms-proxy-init --valid 192:00 -voms cms
# voms-proxy-info

# To run:
# ./make_file_list.sh dataset_name output_file_name

# User inputs:
# - dataset name
# - output file name

# Output:
# - text file containing file paths for dataset

# User inputs:
DATASET=$1
OUTPUT=$2

echo "DATASET: ${DATASET}"
echo "OUTPUT: ${OUTPUT}"

# Check if DATASET is empty.
if [[ -z "$DATASET" ]]; then
    echo "ERROR: DATASET is empty."
    echo "Please provide a dataset name as the first argument."
    exit 1
fi

# Check if OUTPUT is empty.
if [[ -z "$OUTPUT" ]]; then
    echo "ERROR: OUTPUT is empty."
    echo "Please provide an ouput file name as the second argument."
    exit 1
fi

echo "Running dasgoclient query..."

dasgoclient -query="file dataset=${DATASET}" > ${OUTPUT}

NUM_LINES=$(wc -l < ${OUTPUT})

echo "Number of lines in output file: ${NUM_LINES}"

echo "Done!"

