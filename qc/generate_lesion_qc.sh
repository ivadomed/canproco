#!/bin/bash
#
# This script brings the segmentation of the T2w image onto the STIR/PSIR space and generates lesion QC (using SCT) for
# STIR/PSIR images. More context in: https://github.com/ivadomed/canproco/issues/31
#
# Dependencies (versions):
# - SCT (5.8)
#
# Usage:
# sct_run_batch -script qc/generate_lesion_qc.sh -path-data <DATA> -path-output <DATA>_2023-XX-XX -jobs 24

# Manual segmentations or labels should be located under:
# PATH_DATA/derivatives/labels/SUBJECT/ses-0X/anat/

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

# -------------------------------------------------------------------------
# CONVENIENCE FUNCTIONS
# -------------------------------------------------------------------------

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

generate_lesion_qc(){
  ###
  #  This function checks if a manual lesion segmentation exists, then:
  #    - If it does, copy it locally and use it to generate lesion QC
  #    - If it doesn't, skip the subject
  ###
  local image="$1"
  local file_sc="$2"
  # Update global variable with segmentation file name
  FILELESION="${image}_lesion"
  FILELESIONMANUAL="${PATH_DATA}/derivatives/labels/${SUBJECT}/anat/${FILELESION}-manual.nii.gz"
  echo "Looking for manual lesion segmentation: $FILELESIONMANUAL"
  if [[ -e $FILELESIONMANUAL ]]; then
    echo "Found! Using manual lesion segmentation."
    rsync -avzh $FILELESIONMANUAL ${FILELESION}.nii.gz

    # Make sure the lesion is binary
    sct_maths -i ${FILELESION}.nii.gz -bin 0 -o ${FILELESION}_bin.nii.gz
    # Generate lesion QC
    sct_qc -i ${image}.nii.gz -s ${file_sc}.nii.gz -d ${FILELESION}_bin.nii.gz -p sct_deepseg_lesion -plane sagittal -qc ${PATH_QC} -qc-subject ${SUBJECT}
  else
    echo "Lesion segmentation not found. Skipping lesion analysis."
    echo ${FILELESIONMANUAL} >> ${PATH_LOG}/missing_files.log
  fi
}

# Retrieve input params and other params
SUBJECT=$1

# get starting time:
start=`date +%s`

# -------------------------------------------------------------------------
# SCRIPT STARTS HERE
# -------------------------------------------------------------------------
# Display useful info for the log, such as SCT version, RAM and CPU cores available
sct_check_dependencies -short

# Go to folder where data will be copied and processed
cd $PATH_DATA_PROCESSED

# Copy source images
# Note: we use '/./' in order to include the sub-folder 'ses-0X'
rsync -Ravzh $PATH_DATA/./$SUBJECT .

# Go to subject folder for source images
cd ${SUBJECT}/anat

# Define variables
# We do a substitution '/' --> '_' in case there is a subfolder 'ses-0X/'
file="${SUBJECT//[\/]/_}"

# -------------------------------------------------------------------------
# T2w sag
# -------------------------------------------------------------------------
# Add suffix corresponding to contrast
file_t2w=${file}_T2w
# Check if T2w image exists
if [[ -f ${file_t2w}.nii.gz ]];then

    # Reorient the raw image to RPI and resample to 0.8mm isotropic voxel to have the same orientation and resolution
    # as the template SC mask
    mv ${file_t2w}.nii.gz ${file_t2w}_raw.nii.gz
    sct_image -i ${file_t2w}_raw.nii.gz -setorient RPI -o ${file_t2w}_raw_RPI.nii.gz
    sct_resample -i ${file_t2w}_raw_RPI.nii.gz -mm 0.8x0.8x0.8 -o ${file_t2w}_raw_RPI_r.nii.gz

    # Rename _raw_RPI_r file (to be BIDS compliant)
    mv ${file_t2w}_raw_RPI_r.nii.gz ${file_t2w}.nii.gz

    # Spinal cord segmentation
    segment_if_does_not_exist ${file_t2w} 't2' 'deepseg'

fi

# -------------------------------------------------------------------------
# STIR
# -------------------------------------------------------------------------
# Add suffix corresponding to contrast
file_stir=${file}_STIR
# Check if STIR image exists
if [[ -f ${file_stir}.nii.gz ]];then

    # Rename raw file
    mv ${file_stir}.nii.gz ${file_stir}_raw.nii.gz

    # Reorient to RPI
    sct_image -i ${file_stir}_raw.nii.gz -setorient RPI -o ${file_stir}_raw_RPI.nii.gz

    # Rename _raw_RPI file (to be BIDS compliant)
    mv ${file_stir}_raw_RPI.nii.gz ${file_stir}.nii.gz

    # Bring T2w segmentation into STIR space
    sct_register_multimodal -i ${file_t2w}.nii.gz -d ${file_stir}.nii.gz -identity 1 -x nn
    # Bring T2w segmentation into STIR space to be able to run sct_analyze_lesion on the STIR image
    sct_apply_transfo -i ${file_t2w}_seg.nii.gz -d ${file_stir}.nii.gz -w warp_${file_t2w}2${file_stir}.nii.gz -x linear -o ${file_t2w}_seg_reg.nii.gz

    # Generate lesion QC
    generate_lesion_qc ${file_stir} ${file_t2w}_seg_reg

fi

# -------------------------------------------------------------------------
# PSIR
# -------------------------------------------------------------------------
# Add suffix corresponding to contrast
file_psir=${file}_PSIR
# Check if PSIR image exists
if [[ -f ${file_psir}.nii.gz ]];then

    # Rename raw file
    mv ${file_psir}.nii.gz ${file_psir}_raw.nii.gz

    # Reorient to RPI
    sct_image -i ${file_psir}_raw.nii.gz -setorient RPI -o ${file_psir}_raw_RPI.nii.gz

    # Rename _raw_RPI file (to be BIDS compliant)
    mv ${file_psir}_raw_RPI.nii.gz ${file_psir}.nii.gz

    # Bring T2w segmentation into PSIR space
    sct_register_multimodal -i ${file_t2w}.nii.gz -d ${file_psir}.nii.gz -identity 1 -x nn
    # Bring T2w segmentation into PSIR space to be able to run sct_analyze_lesion on the STIR image
    sct_apply_transfo -i ${file_t2w}_seg.nii.gz -d ${file_psir}.nii.gz -w warp_${file_t2w}2${file_psir}.nii.gz -x linear -o ${file_t2w}_seg_reg.nii.gz

    # Generate lesion QC
    generate_lesion_qc ${file_psir} ${file_t2w}_seg_reg

fi

# -------------------------------------------------------------------------
# END
# -------------------------------------------------------------------------

# Display useful info for the log
end=`date +%s`
runtime=$((end-start))
echo
echo "~~~"
echo "SCT version: `sct_version`"
echo "Ran on:      `uname -nsr`"
echo "Duration:    $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
echo "~~~"