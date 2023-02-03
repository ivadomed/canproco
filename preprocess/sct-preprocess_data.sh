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
    # Note: For T2w images, we use sct_deepseg_sc with 2d kernel. Generally, it works better than sct_propseg and sct_deepseg_sc with 3d kernel.
    segment_if_does_not_exist ${file_t2w} 't2' 'deepseg'

    # Vertebral labeling
    label_if_does_not_exist ${file_t2w} ${file_t2w}_seg

    # Compute average cord CSA between C2 and C3
    sct_process_segmentation -i ${file_t2w}_seg.nii.gz -vert 2:3 -vertfile ${file_t2w}_seg_labeled.nii.gz -o ${PATH_RESULTS}/csa-SC_T2w.csv -append 1
    exit

    # MS lesions segmentation
    # TODO - explore why sct_deepseg_lesion produces different results with manually provided centerline
    lesionseg_if_does_not_exist ${file_t2w} 't2'
    # Get centerline from deppeseg SC seg
    sct_get_centerline -i ${file_t2w}_seg.nii.gz -method fitseg
    lesionseg_if_does_not_exist ${file_t2w} 't2' ${file_t2w}_seg_centerline.nii.gz
    exit
fi

# -------------------------------------------------------------------------
# PSIR
# -------------------------------------------------------------------------
# Add suffix corresponding to contrast
file_psir=${file}_PSIR
# Check if PSIR image exists
if [[ -f ${file_psir}.nii.gz ]];then
    # Make sure the image metadata is a valid JSON object
    if [[ ! -s ${file_psir}.json ]]; then
      echo "{}" >> ${file_psir}.json
    fi

    # Rename raw file
    mv ${file_psir}.nii.gz ${file_psir}_raw.nii.gz

    # Reorient to RPI
    sct_image -i ${file_psir}_raw.nii.gz -setorient RPI -o ${file_psir}_raw_RPI.nii.gz

    # Create a symbolic link to _raw_RPI file (to be BIDS compliant)
    mv ${file_psir}_raw_RPI.nii.gz ${file_psir}.nii.gz

    # Spinal cord segmentation
    # Note: For PSIR images, we use sct_propseg. Generally, it works better than sct_deepseg_sc.
    segment_if_does_not_exist ${file_psir} 't1' 'propseg'
fi

# -------------------------------------------------------------------------
# STIR
# -------------------------------------------------------------------------
# Add suffix corresponding to contrast
file_stir=${file}_STIR
# Check if STIR image exists
if [[ -f ${file_stir}.nii.gz ]];then
    # Make sure the image metadata is a valid JSON object
    if [[ ! -s ${file_stir}.json ]]; then
      echo "{}" >> ${file_stir}.json
    fi

    # Rename raw file
    mv ${file_stir}.nii.gz ${file_stir}_raw.nii.gz

    # Reorient to RPI
    sct_image -i ${file_stir}_raw.nii.gz -setorient RPI -o ${file_stir}_raw_RPI.nii.gz

    # Create a symbolic link to _raw_RPI file (to be BIDS compliant)
    mv ${file_stir}_raw_RPI.nii.gz ${file_stir}.nii.gz

    # Spinal cord segmentation
    # Note: For STIR images, we use sct_propseg. Generally, it works better than sct_deepseg_sc.
    segment_if_does_not_exist ${file_stir} 't2' 'propseg'
fi

# ------------------------------------------------------------------------------
# T2star
# ------------------------------------------------------------------------------
# Add suffix corresponding to contrast
file_t2s="${file}_T2star"
# Check if T2star image exists
if [[ -f ${file_t2s}.nii.gz ]];then
    # Make sure the image metadata is a valid JSON object
    if [[ ! -s ${file_t2s}.json ]]; then
      echo "{}" >> ${file_t2s}.json
    fi

    # Rename raw file
    mv ${file_t2s}.nii.gz ${file_t2s}_raw.nii.gz

    # Reorient to RPI
    sct_image -i ${file_t2s}_raw.nii.gz -setorient RPI -o ${file_t2s}_raw_RPI.nii.gz
    # Compute root-mean square across 4th dimension (if it exists)
    sct_maths -i ${file_t2s}_raw_RPI.nii.gz -rms t -o ${file_t2s}_raw_RPI_rms.nii.gz

    # Create a symbolic link to _raw_RPI_rms file (to be BIDS compliant)
    mv ${file_t2s}_raw_RPI_rms.nii.gz ${file_t2s}.nii.gz

    # Spinal cord segmentation
    # Note: For T2star images, we use sct_deepseg_sc
    segment_if_does_not_exist ${file_t2s} 't2s' 'deepseg'

    # TODO - bring vertebral levels from T2w into T2star - so far we do not have T2w_seg_labeled since we are working
    # on manual T2w SC seg corrections
    # sct_register_multimodal -i label_T1w/template/PAM50_levels.nii.gz -d ${file_t2s}.nii.gz -o PAM50_levels2${file_t2s}.nii.gz -identity 1 -x nn

fi
# # ------------------------------------------------------------------------------
# # MTS
# # ------------------------------------------------------------------------------
# # Go to subject folder for source images
# cd ${SUBJECT}/anat
# # Define variables
# # We do a substitution '/' --> '_' in case there is a subfolder 'ses-0X/'
# file="${SUBJECT//[\/]/_}"

# file_t1w="${file}_acq-T1w_MTS"
# file_mton="${file}_acq-MTon_MTS"
# file_mtoff="${file}_acq-MToff_MTS"
# file_mt="${file}_acq-MT_MTS"

# if [[ -e "${file_t1w}.nii.gz" && -e "${file_mton}.nii.gz" && -e "${file_mtoff}.nii.gz" && -e "${file_mt}.nii.gz" ]]; then
#   # Segment lesion (only if it does not exist)
#   lesionseg_if_does_not_exist $file_t1w "t1" ${CENTERLINE_METHOD}
#   file_lesionseg_t1w=$FILELESION

#   # Create mask
#   sct_create_mask -i ${file_t1w}.nii.gz -p centerline,${file_t1w_seg}.nii.gz -size 35mm -o ${file_t1w}_mask.nii.gz
#   # Crop data for faster processing
#   sct_crop_image -i ${file_t1w}.nii.gz -m ${file_t1w}_mask.nii.gz -o ${file_t1w}_crop.nii.gz
#   # file_t1w="${file_t1w}_crop"

#   # Redefine variable for final SC segmentation mask as path changed
#   file_seg_dil_t1w=${PATH_DATA_PROCESSED}/${SUBJECT}/anat/${file_t1w_seg}_dilate
#   # Make sure the first rater metadata is a valid JSON object
#   if [[ ! -s ${file_gt}.json ]]; then
#     echo "{}" >> ${file_gt}.json
#   fi
#   # Crop the manual seg
#   sct_crop_image -i ${file_gt}.nii.gz -m ${file_seg_dil_t1w}.nii.gz -o ${file_gt}_crop.nii.gz
#   # Go back to the root output path
#   cd $PATH_OUTPUT
#   # Create and populate clean data processed folder for training
#   PATH_DATA_PROCESSED_CLEAN="${PATH_DATA_PROCESSED}_clean"
#   # Copy over required BIDs files
#   mkdir -p $PATH_DATA_PROCESSED_CLEAN $PATH_DATA_PROCESSED_CLEAN/${SUBJECT} $PATH_DATA_PROCESSED_CLEAN/${SUBJECT}/anat
#   rsync -avzh $PATH_DATA_PROCESSED/dataset_description.json $PATH_DATA_PROCESSED_CLEAN/
#   rsync -avzh $PATH_DATA_PROCESSED/participants.* $PATH_DATA_PROCESSED_CLEAN/
#   rsync -avzh $PATH_DATA_PROCESSED/README $PATH_DATA_PROCESSED_CLEAN/
#   rsync -avzh $PATH_DATA_PROCESSED/dataset_description.json $PATH_DATA_PROCESSED_CLEAN/derivatives/
#   # For lesion segmentation task, copy SC crops as inputs and lesion annotations as targets
#   rsync -avzh $PATH_DATA_PROCESSED/${SUBJECT}/anat/${file}_crop.nii.gz $PATH_DATA_PROCESSED_CLEAN/${SUBJECT}/anat/${file}.nii.gzROCESSED_CLEAN/${SUBJECT}/anat/${file}.json
#   mkdir -p $PATH_DATA_PROCESSED_CLEAN/derivatives $PATH_DATA_PROCESSED_CLEAN/derivatives/labels $PATH_DATA_PROCESSED_CLEAN/derivatives/labels/${SUBJECT} $PATH_DATA_PROCESSED_CLEAN/derivatives/labels/${SUBJECT}/anat/
#   rsync -avzh $PATH_DATA_PROCESSED/derivatives/labels/${SUBJECT}/anat/${file_gt}_crop.nii.gz $PATH_DATA_PROCESSED_CLEAN/derivatives/labels/${SUBJECT}/anat/${file_gt}.nii.gz
#   rsync -avzh $PATH_DATA_PROCESSED/derivatives/labels/${SUBJECT}/anat/${file_gt}.json $PATH_DATA_PROCESSED_CLEAN/derivatives/labels/${SUBJECT}/anat/${file_gt}.json


#   # Register PD (mtoff)->T1w
#   # Tips: here we only use rigid transformation because both images have very similar sequence parameters. We don't want to use SyN/BSplineSyN to avoid introducing spurious deformations.
#   sct_register_multimodal -i ${file_mtoff}.nii.gz -d ${file_t1w}.nii.gz -dseg ${file_t1w_seg}.nii.gz -param step=1,type=im,algo=rigid,slicewise=1,metric=CC -x spline -qc ${PATH_QC} -qc-subject ${SUBJECT}
#   file_mtoff="${file_mtoff}_reg"
#   # Register MTon->T1w
#   sct_register_multimodal -i ${file_mton}.nii.gz -d ${file_t1w}.nii.gz -dseg ${file_t1w_seg}.nii.gz -param step=1,type=im,algo=rigid,slicewise=1,metric=CC -x spline -qc ${PATH_QC} -qc-subject ${SUBJECT}
#   file_mton="${file_mton}_reg"
#     # Register MT->T1w
#   sct_register_multimodal -i ${file_mt}.nii.gz -d ${file_t1w}.nii.gz -dseg ${file_t1w_seg}.nii.gz -param step=1,type=im,algo=rigid,slicewise=1,metric=CC -x spline -qc ${PATH_QC} -qc-subject ${SUBJECT}
#   file_mton="${file_mton}_reg"
#   # Copy json files to match file basename (it will later be used by sct_compute_mtsat)
#   cp ${SUBJECT}_acq-T1w_MTS.json ${file_t1w}.json
#   cp ${SUBJECT}_acq-MToff_MTS.json ${file_mtoff}.json
#   cp ${SUBJECT}_acq-MTon_MTS.json ${file_mton}.json
#   # # Register template->T1w_ax (using template-T1w as initial transformation)
#   # sct_register_multimodal -i $SCT_DIR/data/PAM50/template/PAM50_t1.nii.gz -iseg $SCT_DIR/data/PAM50/template/PAM50_cord.nii.gz -d ${file_t1w}.nii.gz -dseg ${file_t1w_seg}.nii.gz -param step=1,type=seg,algo=slicereg,metric=MeanSquares,smooth=2:step=2,type=im,algo=syn,metric=CC,iter=5,gradStep=0.5 -initwarp warp_template2T1w.nii.gz -initwarpinv warp_T1w2template.nii.gz
#   # # Rename warping field for clarity
#   # mv warp_PAM50_t12${file_t1w}.nii.gz warp_template2axT1w.nii.gz
#   # mv warp_${file_t1w}2PAM50_t1.nii.gz warp_axT1w2template.nii.gz
#   # # Warp template
#   # sct_warp_template -d ${file_t1w}.nii.gz -w warp_template2axT1w.nii.gz -ofolder label_axT1w -qc ${PATH_QC} -qc-subject ${SUBJECT}
#   # # Compute MTR
#   # sct_compute_mtr -mt0 ${file_mtoff}.nii.gz -mt1 ${file_mton}.nii.gz
#   # # Compute MTsat
#   # sct_compute_mtsat -mt ${file_mton}.nii.gz -pd ${file_mtoff}.nii.gz -t1 ${file_t1w}.nii.gz
#   # # Extract MTR, MTsat and T1 in WM between C2 and C5 vertebral levels
#   # sct_extract_metric -i mtr.nii.gz -f label_axT1w/atlas -l 51 -vert 2:5 -vertfile label_axT1w/template/PAM50_levels.nii.gz -o ${PATH_RESULTS}/MTR.csv -append 1
#   # sct_extract_metric -i mtsat.nii.gz -f label_axT1w/atlas -l 51 -vert 2:5 -vertfile label_axT1w/template/PAM50_levels.nii.gz -o ${PATH_RESULTS}/MTsat.csv -append 1
#   # sct_extract_metric -i t1map.nii.gz -f label_axT1w/atlas -l 51 -vert 2:5 -vertfile label_axT1w/template/PAM50_levels.nii.gz -o ${PATH_RESULTS}/T1.csv -append 1
#   # # Compute MTR in LCST between C2 and C5 vertebral levels
#   # sct_extract_metric -i mtr.nii.gz -f label_axT1w/atlas -l 2,17 -vert 2:5 -vertfile label_axT1w/template/PAM50_levels.nii.gz -o ${PATH_RESULTS}/MTR_LCST.csv -append 1 -combine 1
#   # # Compute MTR in LCST between C2 and C5 vertebral levels
#   # sct_extract_metric -i mtr.nii.gz -f label_axT1w/atlas -l 0,1,15,16 -vert 2:5 -vertfile label_axT1w/template/PAM50_levels.nii.gz -o ${PATH_RESULTS}/MTR_DC.csv -append 1 -combine 1
# else
#   echo "WARNING: MTS dataset is incomplete."
# fi

# Display useful info for the log
end=`date +%s`
runtime=$((end-start))
echo
echo "~~~"
echo "SCT version: `sct_version`"
echo "Ran on:      `uname -nsr`"
echo "Duration:    $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
echo "~~~"







# # # Add suffix corresponding to contrast
# # if [[ ! -f ${file}MT_MTS.nii.gz ]]; then 
# #   file=${file}_acq-MT_MTS
# # elif [[ ! -f ${file}_STIR.nii.gz ]]; then
# #   file=${file}_STIR
# # elif [[ ! -f ${file}_T2star.nii.gz ]]; then
# # file=${file}_T2star
# # elif [[ ! -f ${file}_T2w.nii.gz ]]; 
# # file=${file}_T2w
# # elif [[ ! -f ${file}_PSIR.nii.gz ]]; then
# #   file=${file}_PSIR
# # elif [[ ! -f ${file}T1w_MTS.nii.gz ]]; 
# #   file=${file}_T1w_MTS
# # elif [[ ! -f ${file}T1w.nii.gz ]]; 
# #   file=${file}_T1w
# # elif [[ ! -f ${file}MTon_MTS.nii.gz ]]; then
# #   file=${file}_acq-MTon_MTS
# # else [[ ! -f ${file}MToff_MTS.nii.gz ]];
# #   file=${file}_acq-MToff_MTS
# # fi

# # Make sure the image metadata is a valid JSON object
# if [[ ! -s ${file}.json ]]; then
#   echo "{}" >> ${file}.json
# fi

# # Spinal cord segmentation using the appropriate image contrast
# # for fil in file;do
# #   if [[ "$fil" == *"$T2w"* ]];then
# #     segment_if_does_not_exist ${file} t2 ${CENTERLINE_METHOD}
# #     file_seg="${FILESEG}"
# #   elif [[ "$fil" == *"$T2star"* ]];then
# #     segment_if_does_not_exist ${file} t2s ${CENTERLINE_METHOD}
# #     file_seg="${FILESEG}"
# #   elif [[ "$fil" == *"$T1w"* ]];then
# #     segment_if_does_not_exist ${file} t1 ${CENTERLINE_METHOD}
# #     file_seg="${FILESEG}"
# #   else 
# #     echo "contrast not supported. Skipping"
# #   fi

# segment_if_does_not_exist ${file} t2 ${CENTERLINE_METHOD}
# file_seg="${FILESEG}"

# # get image dims from nifti header info
# ## get header
# ## | = from that output, apply next function / call
# ### -f1 = select 1st col; -d ',' = using ',' as the delimiter 
# ####  grep dim = get all lines that have word 'dim' 
# ##### sed -n '1p' = select first line. See https://riptutorial.com/sed/example/13752/specific-range-of-lines 
# ###### -f2 = select 2nd col; -d '[' = delimited by '['
# ####### $( ) = assign to shell variable
# # dims=$(sct_image -i ${file}.nii.gz -header | cut -f1 -d ',' | grep dim | sed -n '1p' | cut -f2 -d '[')