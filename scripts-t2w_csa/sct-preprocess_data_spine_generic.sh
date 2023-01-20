#!/bin/bash
#
# Compute cross-sectional area (CSA) at C2-C3 vertebral levels from T2w images for spine-generic dataset
#
# Dependencies (versions):
# - SCT (5.8)
#

# Usage:
#     sct_run_batch -c <PATH_TO_REPO>/etc/config_preprocess_data.json
# To run on one sub only; specify this sub in the config_preprocess_data_include.json config file
#     sct_run_batch -c <PATH_TO_REPO>/etc/config_preprocess_data_include.json
# To run on several subs only; specify these subs in the config_preprocess_data_include_list.json config file
#     sct_run_batch -c <PATH_TO_REPO>/etc/config_preprocess_data_include_list.json


# Manual segmentations or labels should be located under:
# PATH_DATA/derivatives/labels/SUBJECT/ses-0X/anat/

# The following global variables are retrieved from the caller sct_run_batch
# but could be overwritten by uncommenting the lines below:
# PATH_DATA_PROCESSED="~/data_processed"
# PATH_RESULTS="~/results"
# PATH_LOG="~/log"
# PATH_QC="~/qc"
#
# Authors: Jan Valosek, Sandrine Bédard, Kiri Stern, Julien Cohen-Adad
# Inspired by the spine-generic process_data.sh script
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

# CONVENIENCE FUNCTIONS
# ======================================================================================================================

segment_if_does_not_exist() {
  ###
  #  This function checks if a manual spinal cord segmentation file already exists, then:
  #    - If it does, copy it locally.
  #    - If it doesn't, perform automatic spinal cord segmentation
  #  This allows you to add manual segmentations on a subject-by-subject basis without disrupting the pipeline.
  ###
  local file="$1"
  local contrast="$2"
  local segmentation_method="$3"  # deepseg or propseg
  # Update global variable with segmentation file name
  FILESEG="${file}_seg"
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
    if [[ $segmentation_method == 'deepseg' ]];then
        sct_deepseg_sc -i ${file}.nii.gz -c ${contrast} -qc ${PATH_QC} -qc-subject ${SUBJECT}
    elif [[ $segmentation_method == 'propseg' ]]; then
        sct_propseg -i ${file}.nii.gz -c ${contrast} -qc ${PATH_QC} -qc-subject ${SUBJECT}
    fi
  fi
}

label_if_does_not_exist(){
  ###
  #  This function checks if a manual labels exists, then:
  #    - If it does, copy it locally and use them to initialize vertebral labeling
  #    - If it doesn't, perform automatic vertebral labeling
  ###
  local file="$1"
  local file_seg="$2"
  # Update global variable with segmentation file name
  FILELABEL="${file}_labels-disc"
  FILELABELMANUAL="${PATH_DATA}/derivatives/labels/${SUBJECT}/anat/${FILELABEL}-manual.nii.gz"
  echo "Looking for manual label: $FILELABELMANUAL"
  if [[ -e $FILELABELMANUAL ]]; then
    echo "Found! Using manual labels."
    rsync -avzh $FILELABELMANUAL ${FILELABEL}.nii.gz
    # Generate labeled segmentation from manual disc labels
    sct_label_vertebrae -i ${file}.nii.gz -s ${file_seg}.nii.gz -discfile ${FILELABEL}.nii.gz -c t2 -qc ${PATH_QC} -qc-subject ${SUBJECT}
  else
    echo "Not found. Proceeding with automatic labeling."
    # Generate vertebral labeling
    sct_label_vertebrae -i ${file}.nii.gz -s ${file_seg}.nii.gz -c t2 -qc ${PATH_QC} -qc-subject ${SUBJECT}
  fi
}

# Retrieve input params and other params
SUBJECT=$1

# get starting time:
start=`date +%s`


# SCRIPT STARTS HERE
# ==============================================================================
# Display useful info for the log, such as SCT version, RAM and CPU cores available
sct_check_dependencies -short

# Go to folder where data will be copied and processed
cd $PATH_DATA_PROCESSED

# Copy BIDS-required files to processed data folder (e.g. list of participants)
if [[ ! -f "participants.tsv" ]]; then
  rsync -avzh $PATH_DATA/participants.tsv .
fi
# Copy list of participants in results folder 
if [[ ! -f "participants.json" ]]; then
  rsync -avzh $PATH_DATA/participants.json .
fi
if [[ ! -f "dataset_description.json" ]]; then
  rsync -avzh $PATH_DATA/dataset_description.json .
fi
#if [[ ! -f "README" ]]; then
#  rsync -avzh $PATH_DATA/README .
#fi

# Copy source images
# Note: we use '/./' in order to include the sub-folder 'ses-0X'
rsync -Ravzh $PATH_DATA/./$SUBJECT/anat/${SUBJECT}_T2w.* .

# # Copy segmentation ground truths (GT)
# mkdir -p derivatives/labels
# rsync -Ravzh $PATH_DATA/derivatives/labels/./$SUBJECT derivatives/labels/.

# Go to subject folder for source images
cd ${SUBJECT}/anat

# Define variables
# We do a substitution '/' --> '_' in case there is a subfolder 'ses-0X/'
file="${SUBJECT//[\/]/_}"

# -------------------------------------------------------------------------
# T2w
# -------------------------------------------------------------------------
# Add suffix corresponding to contrast
file_t2w=${file}_T2w
# Check if T2w image exists
if [[ -f ${file_t2w}.nii.gz ]];then
    # Make sure the image metadata is a valid JSON object
    if [[ ! -s ${file_t2w}.json ]]; then
      echo "{}" >> ${file_t2w}.json
    fi

    # Rename raw file
    mv ${file_t2w}.nii.gz ${file_t2w}_raw.nii.gz

    # Reorient to RPI and resample to 0.8mm isotropic voxel (supposed to be the effective resolution)
    sct_image -i ${file_t2w}_raw.nii.gz -setorient RPI -o ${file_t2w}_raw_RPI.nii.gz
    sct_resample -i ${file_t2w}_raw_RPI.nii.gz -mm 0.8x0.8x0.8 -o ${file_t2w}_raw_RPI_r.nii.gz

    # Rename _raw_RPI_r file (to be BIDS compliant)
    mv ${file_t2w}_raw_RPI_r.nii.gz ${file_t2w}.nii.gz

    # Spinal cord segmentation
    # Note: For T2w images, we use sct_deepseg_sc with 2d kernel. Generally, it works better than sct_propseg and sct_deepseg_sc with 3d kernel.
    segment_if_does_not_exist ${file_t2w} 't2' 'deepseg'

    # Vertebral labeling
    label_if_does_not_exist ${file_t2w} ${file_t2w}_seg

    # Compute average cord CSA between C2 and C3
    sct_process_segmentation -i ${file_t2w}_seg.nii.gz -vert 2:3 -vertfile ${file_t2w}_seg_labeled.nii.gz -o ${PATH_RESULTS}/csa-SC_T2w.csv -append 1
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
