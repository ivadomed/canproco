#!/bin/bash
#
# Register PSIR/STIR contrasts to PAM50 space using the spinal cord segmentation and disc labels
# Then, use the transformation to bring the GT lesion to PAM50 space
#
# NOTE: the commands below assume that GT SC and lesion segmentation and manual disc labels are located under
# derivatives/labels
#
# Usage:
#     sct_run_batch -config config.json
#
#
# Example of config.json:
# {
#  "path_data"   : "<PATH_TO_DATASET>/canproco",
#  "path_output" : "<PATH_TO_DATASET>/canproco_register_to_PAM50_2023-10-21",
#  "script"      : "<PATH_TO_REPO>/canproco/segment_sc_contrast-agnostic/02_register_to_pam50.sh",
#  "jobs"        : 16,
#  "exclude_list"     : "sub-mon118 sub-mon006 ..."
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

# Copy manually created disc labels from derivatives/labels
copy_disc_labels(){
  local file="$1"
  # Update global variable with disc labels file name
  FILELABEL="${file}_labels-disc"
  FILELABELMANUAL="${PATH_DATA}/derivatives/labels/${SUBJECT}/anat/${FILELABEL}.nii.gz"
  echo "Looking for manual disc labels: $FILELABELMANUAL"
  if [[ -e $FILELABELMANUAL ]]; then
    echo "Found! Copying manual disc labels."
    rsync -avzh $FILELABELMANUAL ${FILELABEL}.nii.gz
  else
      echo "File ${FILESEGMANUAL}.nii.gz does not exist" >> ${PATH_LOG}/missing_files.log
      echo "ERROR: Manual disc labels ${FILESEGMANUAL}.nii.gz does not exist. Exiting."
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

    # Create a variable with fname corresponding to the files under derivatives/labels (because we modify the variable
    # for PSIR images)
    file_gt=${file}

    if [[ ${file} =~ "PSIR" ]]; then
      # For PSIR, swap contrast from light cord and dark CSF to dark cord and light CSF to use -c t2 during template registration
      # Context: https://github.com/ivadomed/canproco/issues/46#issuecomment-1752142304
      # TODO: this line will be deleted once the SCT script will include the flag for contrast swapping
      sct_maths -i ${file}.nii.gz -mul -1 -o ${file}_mul.nii.gz
      file=${file}_mul
    fi

    # NOTE: the commands below assume that GT SC, GT lesion segmentation and manual disc labels are located under
    # derivatives/labels
    # Copy GT SC seg
    copy_gt "${file_gt}" "seg"
    file_seg="${file_gt}_seg-manual"

    # Copy GT lesion seg
    copy_gt "${file_gt}" "lesion"
    file_lesion="${file_gt}_lesion-manual"

    # Copy manually created disc labels from derivatives/labels
    copy_disc_labels "${file_gt}"
    file_disc="${file_gt}_labels-disc"

    # Binarize GT lesion segmentation
    sct_maths -i ${file_lesion}.nii.gz -bin 0 -o ${file_lesion}_bin.nii.gz

    # Register to template
    # The flag '-param' is inspired by: Eden, D. et al. Brain, 2019. https://doi.org/10.1093/brain/awy352
    # https://github.com/neuropoly/lesion-mapping/blob/497584281526bbcee2f56bfa9ec076982c7f1503/spinalcord/1_register_data.py#L22
    # Note: we use '-ldisc' to use all disc labels
    sct_register_to_template -i ${file}.nii.gz -s ${file_seg}.nii.gz -ldisc ${file_disc}.nii.gz -c t2 -param step=1,type=seg,algo=centermass,metric=MeanSquares,slicewise=1:step=2,type=seg,algo=bsplinesyn,metric=MeanSquares,slicewise=1,iter=3 -qc ${PATH_QC} -qc-subject ${SUBJECT}

    # Generate additional QC to check the registration
    # Native image in PAM50 space
    # https://github.com/spinalcordtoolbox/spinalcordtoolbox/issues/4166#issuecomment-1773808040
    sct_qc -i ${SCT_DIR}/data/PAM50/template/PAM50_t2.nii.gz -s ${SCT_DIR}/data/PAM50/template/PAM50_cord.nii.gz -d anat2template.nii.gz -p sct_image_stitch -qc ${PATH_QC} -qc-subject ${SUBJECT}
    # Native image in PAM50 space overlaid with PAM50_levels
    # https://github.com/spinalcordtoolbox/spinalcordtoolbox/issues/4166#issuecomment-1773810021
    sct_qc -i anat2template.nii.gz -s ${SCT_DIR}/data/PAM50/template/PAM50_levels.nii.gz -p sct_label_vertebrae -qc ${PATH_QC} -qc-subject ${SUBJECT}

    # Bring GT lesion segmentation to template space
    sct_apply_transfo -i ${file_lesion}_bin.nii.gz -d ${SCT_DIR}/data/PAM50/template/PAM50_t2.nii.gz -w warp_anat2template.nii.gz -x nn -o ${file_lesion}_bin_reg.nii.gz
    # Generate QC (native image in PAM50 space overlaid with GT lesion segmentation in PAM50 space)
    sct_qc -i anat2template.nii.gz -d ${file_lesion}_bin_reg.nii.gz -s ${SCT_DIR}/data/PAM50/template/PAM50_cord.nii.gz -p sct_deepseg_lesion -plane sagittal -qc ${PATH_QC} -qc-subject ${SUBJECT}
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
