#!/bin/bash
#
# Preprocess data.
#
# Dependencies (versions):
# - SCT (5.8)
#
# reference from https://github.com/spine-generic/spine-generic/blob/3aa69465063d823088b2c819d3ad9b81725b2fc8/process_data.sh#L246

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

lesionseg_if_does_not_exist() {
 ###
 #  This function checks if a manual lesions segmentation file already exists, then:
 #    - If it does, copy it locally.
 #    - If it doesn't, perform automatic lesions segmentation
 #  This allows you to add manual segmentations on a subject-by-subject basis without disrupting the pipeline.
 ###
 local file="$1"
 local contrast="$2"
 local centerline="$3"
 # Update global variable with segmentation file name
 FILELESION="${file}_lesionseg"
 FILELESIONMANUAL="${PATH_DATA}/derivatives/labels/${SUBJECT}/anat/${FILELESION}-manual.nii.gz"
 echo
 echo "Looking for manual segmentation: $FILELESIONMANUAL"
 if [[ -e $FILELESIONMANUAL ]]; then
   echo "Found! Using manual segmentation."
   rsync -avzh $FILELESIONMANUAL ${FILELESION}.nii.gz
   sct_qc -i ${file}.nii.gz -s ${FILELESION}.nii.gz -p sct_deepseg_lesion -qc ${PATH_QC} -qc-subject ${SUBJECT}
 else
   echo "Not found. Proceeding with automatic segmentation."
   # Lesions segmentation
   # Note, sct_deepseg_lesion does not have QC implemented yet, see: https://github.com/spinalcordtoolbox/spinalcordtoolbox/issues/3803
   if [[ ${centerline} == "" ]];then
      sct_deepseg_lesion -i ${file}.nii.gz -c ${contrast} -ofolder without_centerline
    else
      sct_deepseg_lesion -i ${file}.nii.gz -centerline file -file_centerline ${centerline} -c ${contrast} -ofolder with_centerline
    fi

 fi
}

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
  FILELABEL="${file}_labels"
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

# Copy source images
# Note: we use '/./' in order to include the sub-folder 'ses-0X'
rsync -Ravzh $PATH_DATA/./$SUBJECT .

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
    segment_if_does_not_exist ${file_t2w} 't2' 'deepseg'
fi

# -------------------------------------------------------------------------
# Co-register other contrast to T2w
# -------------------------------------------------------------------------
# Initialize filenames
file_stir="${file}_STIR"
file_psir="${file}_PSIR"
file_t2s="${file}_T2star"
file_t1_mts="${file}_acq-T1w_MTS"
file_mton_mts="${file}_acq-MTon_MTS"
file_mtoff_mts="${file}_acq-MToff_MTS"

contrasts=($file_stir $file_psir $file_t2s $file_t1_mts $file_mton_mts $file_mtoff_mts)


# Prepare T2star images
# Check if T2star image exists
if [[ -f ${file_t2s}.nii.gz ]];then
    # Rename raw file
    mv ${file_t2s}.nii.gz ${file_t2s}_raw.nii.gz
    file_t2s="${file_t2s}_raw"
    # Compute root-mean square across 4th dimension (if it exists), corresponding to all echoes in Philips scans.
    sct_maths -i ${file_t2s}.nii.gz -rms t -o ${file_t2s}_rms.nii.gz
    file_t2s="${file_t2s}_rms"
    # Rename _rms file
    mv ${file_t2s}.nii.gz ${file}_T2star.nii.gz

    # Spinal cord segmentation
    # Note: For T2star images, we use sct_deepseg_sc
    segment_if_does_not_exist ${file}_T2star 't2s' 'deepseg'
    # TODO - bring vertebral levels from T2w into T2star
fi


# Loop across contrasts
for contrast in "${contrasts[@]}"; do
    # Check if contrast exists
    if [[ -f ${contrast}.nii.gz ]];then
        # Bring contrast to T2w space
        sct_register_multimodal -i ${contrast}.nii.gz -d ${file_t2w}.nii.gz -o ${contrast}2${file_t2w}.nii.gz -identity 1 -x nn

        FILESEG_T2w=${file_t2w}_seg
        # Create QC report to assess registration quality
        # Note: registration quality is assessed by comparing the ${contrast} image to the T2w SC seg
        sct_qc -i ${contrast}2${file_t2w}.nii.gz -s ${FILESEG_T2w}.nii.gz -p sct_get_centerline -qc ${PATH_QC} -qc-subject ${SUBJECT}
        sct_qc -i ${contrast}2${file_t2w}.nii.gz -s ${FILESEG_T2w}.nii.gz -p sct_label_vertebrae -qc ${PATH_QC} -qc-subject ${SUBJECT}
   fi
done


# -------------------------------------------------------------------------
# PSIR
# -------------------------------------------------------------------------
file_psir="${file}_PSIR"
# Check if PSIR image exists
if [[ -f ${file_psir}.nii.gz ]];then
    # Spinal cord segmentation
    # Try sct_propseg
    sct_propseg -i ${file_psir}.nii.gz -c t1 -o ${file_psir}_seg_propseg_t1.nii.gz -qc ${PATH_QC} -qc-subject ${SUBJECT}
    # Try sct_deepseg_sc 2d vs 3d kernels
    sct_deepseg_sc -i ${file_psir}.nii.gz -c t1 -o ${file_psir}_seg_deepseg_sc_t1.nii.gz -qc ${PATH_QC} -qc-subject ${SUBJECT}
    sct_deepseg_sc -i ${file_psir}.nii.gz -c t1 -kernel 3d -o ${file_psir}_seg_deepseg_sc_t1_kernel3d.nii.gz -qc ${PATH_QC} -qc-subject ${SUBJECT}
fi

# -------------------------------------------------------------------------
# STIR
# -------------------------------------------------------------------------
file_stir="${file}_STIR"
# Check if STIR image exists
if [[ -f ${file_stir}.nii.gz ]];then
    # Spinal cord segmentation
    # Note: we use the T2w contrast for STIR segmentation
    segment_if_does_not_exist $file_stir 't2' 'deepseg'
fi

# # ------------------------------------------------------------------------------
# # MT
# # ------------------------------------------------------------------------------
if [[ -f ${file_t1_mts}.nii.gz ]];then
    # Spinal cord segmentation
    segment_if_does_not_exist ${file_t1_mts} 't1' 'deepseg'
fi

if [[ -f ${file_mton_mts}.nii.gz ]];then
    # Spinal cord segmentation
    segment_if_does_not_exist ${file_mton_mts} 't2s' 'deepseg'
fi

if [[ -f ${file_mtoff_mts}.nii.gz ]];then
    # Spinal cord segmentation
    segment_if_does_not_exist ${file_mtoff_mts} 't1' 'deepseg'
fi

# # ------------------------------------------------------------------------------
# Display useful info for the log
end=`date +%s`
runtime=$((end-start))
echo
echo "~~~"
echo "SCT version: `sct_version`"
echo "Ran on:      `uname -nsr`"
echo "Duration:    $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
echo "~~~"
