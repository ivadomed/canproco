#!/bin/bash
#
# Segment SC using the contrast-agnostic MONAI model from PSIR/STIR contrast and perform vertebral labeling
#
# Usage:
#     sct_run_batch -config config.json
#
# Note: conda environment with MONAI is required to run this script:
#     conda create -n monai python=3.9
#     conda activate monai
#     pip install -r segment_sc_contrast-agnostic/requirements.txt
#
# Example of config.json:
# {
#  "path_data"   : "<PATH_TO_DATASET>/canproco",
#  "path_output" : "<PATH_TO_DATASET>/canproco_contrast-agnostic_2023-10-06",
#  "script"      : "<PATH_TO_REPO>/canproco/segment_sc_contrast-agnostic/segment_sc_contrast-agnostic.sh",
#  "jobs"        : 16,
#  "exclude"     : "sub-mon118_ses-M0",
#  "script_args" : "<PATH_TO_REPO>/segment_sc_contrast-agnostic/run_inference_single_image.py <PATH_TO_CONTRAST_AGNOSTIC_MODEL>"
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

  # Run SC segmentation
  # TODO: the following call will be replaced by the SCT script
  python ${PATH_SCRIPT} --path-img ${file}.nii.gz --path-out ./ --chkp-path ${PATH_MODEL}

  # Generate sagittal QC report (https://github.com/ivadomed/canproco/issues/37#issuecomment-1644497220)
  sct_qc -i ${file}.nii.gz -s ${file}_pred.nii.gz -d ${file}_pred.nii.gz -p sct_deepseg_lesion -plane sagittal -qc ${PATH_QC} -qc-subject ${SUBJECT}
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

# Copy GT lesion segmentation
copy_gt(){
  local file="$1"
  local type="$2"     # seg or lesion
  # Construct file name to GT lesion segmentation located under derivatives/labels
  FILESEGMANUAL="${PATH_DATA}/derivatives/labels/${SUBJECT}/anat/${file}_${type}-manual.nii.gz"
  echo ""
  echo "Looking for manual segmentation: $FILESEGMANUAL"
  if [[ -e $FILESEGMANUAL ]]; then
      echo "Found! Copying ..."
      rsync -avzh $FILESEGMANUAL ${file}_${type}-manual.nii.gz
  else
      echo "File ${FILESEGMANUAL}.nii.gz does not exist" >> ${PATH_LOG}/missing_files.log
      echo "ERROR: Manual GT segmentation ${FILESEGMANUAL}.nii.gz does not exist. Exiting."
      exit 1
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

    # Copy GT lesion seg
    # Note: we do copy here before we modify the file variable for PSIR
    copy_gt "${file}" "lesion"

    if [[ ${file} =~ "PSIR" ]]; then
      # For PSIR, swap contrast from light cord and dark CSF to dark cord and light CSF
      # Context: https://github.com/ivadomed/canproco/issues/46#issuecomment-1752142304
      # TODO: this line will be deleted once the SCT script will include the flag for contrast swapping
      sct_maths -i ${file}.nii.gz -mul -1 -o ${file}_mul.nii.gz
      file=${file}_mul
    fi

    # Segment SC using the contrast agnostic MONAI model
    segment_sc_monai "${file}"

    # Binarize GT lesion segmentation (sct_analyze_lesion requires binary mask)
    sct_maths -i ${file}_lesion-manual.nii.gz -bin 0 -o ${file}_lesion-manual_bin.nii.gz
    # Analyze GT MS lesion
    sct_analyze_lesion -m ${file}_lesion-manual_bin.nii.gz -s ${file}_seg_monai.nii.gz -ofolder ${PATH_RESULTS}

    # Perform vertebral labeling
    # STIR and PSIR_mul (cord dark; CSF bright) --> T2w
    label_if_does_not_exist "${file}" "${file}_seg_monai" "t2"

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
