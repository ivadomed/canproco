#!/bin/bash
#
# Process spine-generic data to compute T2w C2-C3 CSA
# Note: since T2w manual disc labels do not contain labels 1 and 2, PAM50 labels from T1w is used. It means that
# also T1w images have to be processed.
#
# Dependencies:
# - SCT
#
# Inspired by https://github.com/spine-generic/spine-generic/blob/master/process_data.sh

# Usage:
#     sct_run_batch -c <config_file>.json

# Manual segmentations or labels should be located under:
# PATH_DATA/derivatives/labels/SUBJECT/anat/

# The following global variables are retrieved from the caller sct_run_batch
# but could be overwritten by uncommenting the lines below:
# PATH_DATA_PROCESSED="~/data_processed"
# PATH_RESULTS="~/results"
# PATH_LOG="~/log"
# PATH_QC="~/qc"
#
# Authors: Jan Valosek
#

# Uncomment for full verbose
set -x

# Immediately exit if error
set -e -o pipefail

# Exit if user presses CTRL+C (Linux) or CMD+C (OSX)
trap "echo Caught Keyboard Interrupt within script. Exiting now.; exit" INT

# Print retrieved variables from sct_run_batch to the log (to allow easier debug)
echo "Retrieved variables from from the caller sct_run_batch:"
echo "PATH_DATA: ${PATH_DATA}"
echo "PATH_DATA_PROCESSED: ${PATH_DATA_PROCESSED}"
echo "PATH_RESULTS: ${PATH_RESULTS}"
echo "PATH_LOG: ${PATH_LOG}"
echo "PATH_QC: ${PATH_QC}"

# Retrieve input params
SUBJECT=$1

# get starting time:
start=`date +%s`

# CONVENIENCE FUNCTIONS
# Check if manual label already exists. If it does, copy it locally. If it does
# not, perform labeling.
label_if_does_not_exist(){
  local file="$1"
  local file_seg="$2"
  # Update global variable with segmentation file name
  FILELABEL="${file/_RPI_r/}_labels-disc"
  FILELABELMANUAL="${PATH_DATA}/derivatives/labels/${SUBJECT}/anat/${FILELABEL}-manual.nii.gz"
  echo "Looking for manual label: $FILELABELMANUAL"
  if [[ -e $FILELABELMANUAL ]]; then
    echo "Found! Using manual labels."
    rsync -avzh $FILELABELMANUAL ${FILELABEL}.nii.gz
  else
    echo "Not found. Proceeding with automatic labeling."
    echo "Missing labels ${file}" >> ${PATH_LOG}/missing_labels.txt
    exit
  fi
}

# Check if manual segmentation already exists. If it does, copy it locally. If
# it does not, perform seg.
segment_if_does_not_exist(){
  local file="$1"
  local contrast="$2"
  # Update global variable with segmentation file name
  FILESEG="${file/_RPI_r/}_seg"
  FILESEGMANUAL="${PATH_DATA}/derivatives/labels/${SUBJECT}/anat/${FILESEG}-manual.nii.gz"
  echo
  echo "Looking for manual segmentation: $FILESEGMANUAL"
  if [[ -e $FILESEGMANUAL ]]; then
    echo "Found! Using manual segmentation."
    rsync -avzh $FILESEGMANUAL ${FILESEG}.nii.gz
    sct_qc -i ${file}.nii.gz -s ${FILESEG}.nii.gz -p sct_deepseg_sc -qc ${PATH_QC} -qc-subject ${SUBJECT}
  else
    echo "Not found. Proceeding with automatic segmentation."
    # Segment spinal cord
    sct_deepseg_sc -i ${file}.nii.gz -c $contrast -qc ${PATH_QC} -qc-subject ${SUBJECT}
  fi
}

# SCRIPT STARTS HERE
# ==============================================================================
# Display useful info for the log, such as SCT version, RAM and CPU cores available
sct_check_dependencies -short

# Go to folder where data will be copied and processed
cd $PATH_DATA_PROCESSED
# Copy list of participants in processed data folder
if [[ ! -f "participants.tsv" ]]; then
  rsync -avzh $PATH_DATA/participants.tsv .
fi
# Copy list of participants in results folder (used by spine-generic scripts)
if [[ ! -f $PATH_RESULTS/"participants.tsv" ]]; then
  rsync -avzh $PATH_DATA/participants.tsv $PATH_RESULTS/"participants.tsv"
fi

# Copy source T1w and T2w images
# Note: we use '/./' in order to include the sub-folder 'ses-0X'
rsync -Ravzh ${PATH_DATA}/./${SUBJECT}/anat/${SUBJECT}_T[12]w.* .
# Go to anat folder where all structural data are located
cd ${SUBJECT}/anat/

# T1w
# ------------------------------------------------------------------------------
file_t1="${SUBJECT}_T1w"
# Reorient to RPI and resample to 1mm iso (supposed to be the effective resolution)
sct_image -i ${file_t1}.nii.gz -setorient RPI -o ${file_t1}_RPI.nii.gz
sct_resample -i ${file_t1}_RPI.nii.gz -mm 1x1x1 -o ${file_t1}_RPI_r.nii.gz
file_t1="${file_t1}_RPI_r"
# Segment spinal cord (only if it does not exist)
segment_if_does_not_exist $file_t1 "t1"
file_t1_seg=$FILESEG
# Create mid-vertebral levels in the cord (only if it does not exist)
label_if_does_not_exist ${file_t1} ${file_t1_seg}
file_label=$FILELABEL
# Register to PAM50 template
sct_register_to_template -i ${file_t1}.nii.gz -s ${file_t1_seg}.nii.gz -l ${file_label}.nii.gz -c t1 -param step=1,type=seg,algo=centermassrot:step=2,type=seg,algo=syn,slicewise=1,smooth=0,iter=5:step=3,type=im,algo=syn,slicewise=1,smooth=0,iter=3 -qc ${PATH_QC} -qc-subject ${SUBJECT}
# Rename warping fields for clarity
mv warp_template2anat.nii.gz warp_template2T1w.nii.gz
mv warp_anat2template.nii.gz warp_T1w2template.nii.gz
# Warp template without the white matter atlas (we don't need it at this point)
sct_warp_template -d ${file_t1}.nii.gz -w warp_template2T1w.nii.gz -a 0 -ofolder label_T1w
# Generate QC report to assess vertebral labeling
sct_qc -i ${file_t1}.nii.gz -s label_T1w/template/PAM50_levels.nii.gz -p sct_label_vertebrae -qc ${PATH_QC} -qc-subject ${SUBJECT}

# T2w
# ------------------------------------------------------------------------------
file_t2="${SUBJECT}_T2w"
# Reorient to RPI and resample to 0.8mm iso (supposed to be the effective resolution)
sct_image -i ${file_t2}.nii.gz -setorient RPI -o ${file_t2}_RPI.nii.gz
sct_resample -i ${file_t2}_RPI.nii.gz -mm 0.8x0.8x0.8 -o ${file_t2}_RPI_r.nii.gz
file_t2="${file_t2}_RPI_r"
# Segment spinal cord (only if it does not exist)
segment_if_does_not_exist $file_t2 "t2"
file_t2_seg=$FILESEG

# Bring vertebral level into T2 space
sct_register_multimodal -i label_T1w/template/PAM50_levels.nii.gz -d ${file_t2_seg}.nii.gz -o PAM50_levels2${file_t2}.nii.gz -identity 1 -x nn

# Generate QC report to assess vertebral labeling
sct_qc -i ${file_t2}.nii.gz -s PAM50_levels2${file_t2}.nii.gz -p sct_label_vertebrae -qc ${PATH_QC} -qc-subject ${SUBJECT}

# Compute average cord CSA between C2 and C3
sct_process_segmentation -i ${file_t2_seg}.nii.gz -vert 2:3 -vertfile PAM50_levels2${file_t2}.nii.gz -o ${PATH_RESULTS}/csa-SC_T2w.csv -append 1

# ------------------------------------------------------------------------------

# Display useful info for the log
end=`date +%s`
runtime=$((end-start))
echo
echo "~~~"
echo "SCT version: `sct_version`"
echo "Ran on:      `uname -nsr`"
echo "Duration:    $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
echo "~~~"