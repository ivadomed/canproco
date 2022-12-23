#!/bin/bash
#
# Create manual MS lesion segmentation from PSIR images.
#
# Usage:
#   ./psir_correct_lesions.sh <SUBJECT>
#
# Manual segmentations created by this script will be located under:
# PATH_DATA_PROCESSED/derivatives/labels/SUBJECT/anat/
#
# Authors: Sandrine Bédard, Jan Valosek, Michelle Chen, Julien Cohen-Adad
#
# How to use: sct_run_batch -jobs 1 -path-data <> -path-output <> -script psir_correct_lesions.sh

# TODO: this script could be merged with manual_correction.py

set -x
# Immediately exit if error
set -e -o pipefail

# Exit if user presses CTRL+C (Linux) or CMD+C (OSX)
trap "echo Caught Keyboard Interrupt within script. Exiting now.; exit" INT

# Retrieve input params
SUBJECT=$1

# Save script path
PATH_SCRIPT=$PWD

# get starting time:
start=`date +%s`

# SCRIPT STARTS HERE
# ==============================================================================
# Display useful info for the log, such as SCT version, RAM and CPU cores available
sct_check_dependencies -short

# Go to folder where data will be copied and processed
cd ${PATH_DATA_PROCESSED}
# Copy list of participants in processed data folder
if [[ ! -f "participants.tsv" ]]; then
  rsync -avzh $PATH_DATA/participants.tsv .
fi
# Copy list of participants in resutls folder
if [[ ! -f $PATH_RESULTS/"participants.tsv" ]]; then
  rsync -avzh $PATH_DATA/participants.tsv $PATH_RESULTS/"participants.tsv"
fi
# Copy source images
rsync -Ravzh $PATH_DATA/./$SUBJECT .
# Go to anat folder where all structural data are located
cd ${SUBJECT}/anat/

# PSIR
# ------------------------------------------------------------------------------
# Define variables
# We do a substitution '/' --> '_' in case there is a subfolder 'ses-0X/'
file="${SUBJECT//[\/]/_}_PSIR"

if [[ -e ${file}.nii.gz ]]; then
    echo "Enter your name (Firstname Lastname). It will be used to generate a json sidecar with each corrected file: "
    read -r name_for_json
    # Create a subject folder under the derivatives
    mkdir -p ${PATH_DATA_PROCESSED}/derivatives/labels/${SUBJECT}/anat/
    # create an empty mask
    sct_create_mask -i ${file}.nii.gz -o ${PATH_DATA_PROCESSED}/derivatives/labels/${SUBJECT}/anat/${file}_lesion-manual.nii.gz -p center -size 0
    # convert mask to proper datatype
    sct_image -i ${PATH_DATA_PROCESSED}/derivatives/labels/${SUBJECT}/anat/${file}_lesion-manual.nii.gz -type uint16 -o ${PATH_DATA_PROCESSED}/derivatives/labels/${SUBJECT}/anat/${file}_lesion-manual.nii.gz
    fsleyes ${file}.nii.gz ${PATH_DATA_PROCESSED}/derivatives/labels/${SUBJECT}/anat/${file}_lesion-manual.nii.gz
    # create a json sidecar (4 spaces are used for indentation)
    cur_time=$(python <<< "import time;print(time.strftime('%Y-%m-%d %H:%M:%S'))")
    echo -e "{\n    \"Author\": \"${name_for_json}\",\n    \"Label\": \"lesion-manual\",\n    \"Date\": \"${cur_time}\"\n}" > ${PATH_DATA_PROCESSED}/derivatives/labels/${SUBJECT}/anat/${file}_lesion-manual.json
fi


# Display useful info for the log
end=`date +%s`
runtime=$((end-start))
echo
echo "~~~"
echo "SCT version: `sct_version`"
echo "Ran on:      `uname -nsr`"
echo "Duration:    $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
echo "~~~"
