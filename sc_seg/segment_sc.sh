#!/bin/bash
#
# Segment SC using the contrast-agnostic MONAI model from PSIR contrast and perform vertebral labeling
#
# Usage:
#     sct_run_batch -config config.json
#
# Note: conda environment with MONAI is required to run this script:
#     conda create -n monai python=3.8
#     conda activate monai
#     pip install -r sc_seg/requirements.txt
#
# Example of config.json:
# {
#  "path_data"   : "<PATH_TO_DATASET>/canproco",
#  "path_output" : "<PATH_TO_DATASET>/canproco_contrast-agnostic_2023-10-06",
#  "script"      : "<PATH_TO_REPO>/canproco/sc_seg/segment_sc.sh",
#  "jobs"        : 16,
#  "exclude"     : "sub-mon118_ses-M0",
#  "script_args" : "<PATH_TO_SCRIPT>/run_inference_single_image.py <PATH_TO_CONTRAST_AGNOSTIC_MODEL>"
# }
#
# The following global variables are retrieved from the caller sct_run_batch
# but could be overwritten by uncommenting the lines below:
# PATH_DATA_PROCESSED="~/data_processed"
# PATH_RESULTS="~/results"
# PATH_LOG="~/log"
# PATH_QC="~/qc"
#
# Author: Jan Valosek
#

# Uncomment for full verbose
set -x

# Immediately exit if error
set -e -o pipefail

# Exit if user presses CTRL+C (Linux) or CMD+C (OSX)
trap "echo Caught Keyboard Interrupt within script. Exiting now.; exit" INT

# Print retrieved variables from the sct_run_batch script to the log (to allow easier debug)
echo "Retrieved variables from from the caller sct_run_batch:"
echo "PATH_DATA: ${PATH_DATA}"
echo "PATH_DATA_PROCESSED: ${PATH_DATA_PROCESSED}"
echo "PATH_RESULTS: ${PATH_RESULTS}"
echo "PATH_LOG: ${PATH_LOG}"
echo "PATH_QC: ${PATH_QC}"

SUBJECT=$1
PATH_SCRIPT=$2
PATH_MODEL=$3

echo "SUBJECT: ${SUBJECT}"
echo "PATH_SCRIPT: ${PATH_SCRIPT}"
echo "PATH_MODEL: ${PATH_MODEL}"

# ------------------------------------------------------------------------------
# CONVENIENCE FUNCTIONS
# ------------------------------------------------------------------------------

# Segment spinal cord using contrast agnostic MONAI model
segment_sc_monai(){
  local file="$1"

  FILESEG="${file}_seg_monai"

  # Get the start time
  start_time=$(date +%s)
  # Run SC segmentation
  python ${PATH_SCRIPT} -i ${file}.nii.gz -o ${FILESEG}.nii.gz -path-model ${PATH_MODEL}
  # Get the end time
  end_time=$(date +%s)
  # Calculate the time difference
  execution_time=$(python3 -c "print($end_time - $start_time)")
  echo "${FILESEG},${execution_time}" >> ${PATH_RESULTS}/execution_time.csv

  # Generate QC report
  sct_qc -i ${file}.nii.gz -s ${FILESEG}.nii.gz -p sct_deepseg_sc -qc ${PATH_QC} -qc-subject ${SUBJECT}
}

# Check if manual label already exists. If it does, copy it locally. If it does
# not, perform labeling.
label_if_does_not_exist(){
  local file="$1"
  local file_seg="$2"
  local contrast="$3"
  # Update global variable with segmentation file name
  FILELABEL="${file}_labels"
  FILELABELMANUAL="${PATH_DATA}/derivatives/labels/${SUBJECT}/anat/${FILELABEL}-manual.nii.gz"
  echo "Looking for manual label: $FILELABELMANUAL"
  if [[ -e $FILELABELMANUAL ]]; then
    echo "Found! Using manual labels."
    rsync -avzh $FILELABELMANUAL ${FILELABEL}.nii.gz
  else
    echo "Not found. Proceeding with automatic labeling."
    # Generate labeled segmentation
    sct_label_vertebrae -i ${file}.nii.gz -s ${file_seg}.nii.gz -c ${contrast} -qc ${PATH_QC} -qc-subject ${SUBJECT}
  fi
}

# ------------------------------------------------------------------------------
# SCRIPT STARTS HERE
# ------------------------------------------------------------------------------
# get starting time:
start=`date +%s`

# Display useful info for the log, such as SCT version, RAM and CPU cores available
sct_check_dependencies -short

# Go to folder where data will be copied and processed
cd $PATH_DATA_PROCESSED

# Copy source PSIR/STIR images
# Note: we use '/./' in order to include the sub-folder 'ses-M0'
# We do a substitution '/' --> '_' in case there is a subfolder 'ses-M0/'
rsync -Ravzh ${PATH_DATA}/./${SUBJECT}/anat/${SUBJECT//[\/]/_}_*IR.* .

# Go to subject folder for source images
cd ${SUBJECT}/anat

# ------------------------------------------------------------------------------
# PSIR/STIR
# ------------------------------------------------------------------------------

# We do a substitution '/' --> '_' in case there is a subfolder 'ses-M0/'
# Calgary has STIR
if [[ ${SUBJECT} =~ "cal" ]]; then
  file="${SUBJECT//[\/]/_}"_STIR
# Other sites have PSIR
else
  file="${SUBJECT//[\/]/_}"_PSIR
fi

# Check if file exists
if [[ ! -e ${file}.nii.gz ]]; then
    echo "File ${file}.nii.gz does not exist" >> ${PATH_LOG}/missing_files.log
    echo "ERROR: File ${file}.nii.gz does not exist. Exiting."
    exit 1
else
    # Segment SC using the contrast agnostic MONAI model
    segment_sc_monai "${file}"

    # Perform vertebral labeling
    if [[ ${file} =~ *"PSIR"* ]]; then
      label_if_does_not_exist "${file}" "${file}_seg_monai" "t1"
    else
      label_if_does_not_exist "${file}" "${file}_seg_monai" "t2"
    fi
fi

# ------------------------------------------------------------------------------
# End
# ------------------------------------------------------------------------------

# Display results (to easily compare integrity across SCT versions)
end=`date +%s`
runtime=$((end-start))
echo
echo "~~~"
echo "SCT version: `sct_version`"
echo "Ran on:      `uname -nsr`"
echo "Duration:    $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
echo "~~~"
