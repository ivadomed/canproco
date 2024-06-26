#!/bin/bash
# 
# Process data for the CanProCo project. This script does the following:
# - Compute CSA in whole spinal cord between C2-C4
# - Compute MTR in whole spinal cord, white matter, gray matter, dorsal columns, lateral funiculi, ventral funiculi between C2-C4
# 
# Usage: sh spine.sh <ABSOLUTE_PATH_SUBJECT> <qc>
# Example command line: sh spine.sh /Users/msresearch/Downloads/canproco/data/CAN-01-CON-039-M0 /Users/msresearch/Downloads/canproco/data/qc
# 
# Tested with Spinal Cord Toolbox - Version 5.7
# 
# Author: Lisa Eunyoung Lee

SUBJECT=$1
qc=$2

# Display start time
echo "${SUBJECT}: spine.sh started at $(date +%x_%r)";
start=$(date +%s);

# T2
# ===========================================================================================

cd ${SUBJECT}/t2;

# Segment spinal cord

sct_image -i t2.nii.gz -setorient RPI -o t2.nii.gz;

sct_deepseg_sc -i t2.nii.gz -c t2 -qc $qc;

# Label vertebral levels

sct_label_vertebrae -i t2.nii.gz -s t2_seg.nii.gz -c t2 -qc $qc;

sct_label_utils -i t2_seg_labeled.nii.gz -vert-body 3,5 -o t2_labels_vert.nii.gz;

sct_qc -i t2.nii.gz -s t2_labels_vert.nii.gz -p sct_label_utils -qc $qc;

# Register t2 to template space then, warp template to t2	

sct_register_to_template -i t2.nii.gz -s t2_seg.nii.gz -l t2_labels_vert.nii.gz -c t2 -qc $qc;

sct_warp_template -d t2.nii.gz -w warp_template2anat.nii.gz -a 0 -qc $qc;

# Compute cross-sectional area averaged between C2-C4

sct_process_segmentation -i t2_seg.nii.gz -vert 2:4 -vertfile ./label/template/PAM50_levels.nii.gz -o csa_c2c4.csv -append 0;


# MT
# ===========================================================================================

cd ${SUBJECT}/mt;

# Segment spinal cord

sct_image -i mt_t1.nii.gz -setorient RPI -o mt_t1.nii.gz;

sct_deepseg_sc -i mt_t1.nii.gz -c t1 -qc $qc;

# Create a mask

sct_create_mask -i mt_t1.nii.gz -p centerline,mt_t1_seg.nii.gz -size 35mm -f cylinder -o mt_t1_mask.nii.gz;

# Register template -> mt_t1 then, warp template to mt_t1 

sct_register_multimodal -i $SCT_DIR/data/PAM50/template/PAM50_t1.nii.gz -iseg $SCT_DIR/data/PAM50/template/PAM50_cord.nii.gz -d mt_t1.nii.gz -dseg mt_t1_seg.nii.gz -m mt_t1_mask.nii.gz -param step=1,type=seg,algo=slicereg,metric=MeanSquares,smooth=2:step=2,type=im,algo=syn,metric=CC,iter=5,gradStep=0.5 -initwarp ${SUBJECT}/t2/warp_template2anat.nii.gz -initwarpinv ${SUBJECT}/t2/warp_anat2template.nii.gz -owarp warp_template2mt.nii.gz -qc $qc;

sct_warp_template -d mt_t1.nii.gz -w warp_template2mt.nii.gz -a 1 -qc $qc;

# Register mt0 -> mt_t1

sct_register_multimodal -i mt0.nii.gz -d mt_t1.nii.gz -dseg mt_t1_seg.nii.gz -m mt_t1_mask.nii.gz -param step=1,type=im,algo=slicereg,metric=CC -x spline -o mt0_reg.nii.gz -qc $qc;

# Register mt1 -> mt_t1

sct_register_multimodal -i mt1.nii.gz -d mt_t1.nii.gz -dseg mt_t1_seg.nii.gz -m mt_t1_mask.nii.gz -param step=1,type=im,algo=slicereg,metric=CC -x spline -o mt1_reg.nii.gz -qc $qc;

# Compute MTR

sct_compute_mtr -mt0 mt0_reg.nii.gz -mt1 mt1_reg.nii.gz;

# Extract MTR in whole cord, GM, WM, WM regions (DC, LF, VF, CST) averaged between C2-C4

mkdir -p ${SUBJECT}/mt/csv;

sct_extract_metric -i mtr.nii.gz -f label/atlas -method map -l 50 -vert 2:4 -vertfile label/template/PAM50_levels.nii.gz -z 3:18 -o csv/mtr_cord.csv -append 0;
sct_extract_metric -i mtr.nii.gz -f label/atlas -method map -l 52 -vert 2:4 -vertfile label/template/PAM50_levels.nii.gz -z 3:18 -o csv/mtr_gm.csv -append 0;
sct_extract_metric -i mtr.nii.gz -f label/atlas -method map -l 51 -vert 2:4 -vertfile label/template/PAM50_levels.nii.gz -z 3:18 -o csv/mtr_wm.csv -append 0;
sct_extract_metric -i mtr.nii.gz -f label/atlas -method map -l 53 -vert 2:4 -vertfile label/template/PAM50_levels.nii.gz -z 3:18 -o csv/mtr_dc.csv -append 0;
sct_extract_metric -i mtr.nii.gz -f label/atlas -method map -l 54 -vert 2:4 -vertfile label/template/PAM50_levels.nii.gz -z 3:18 -o csv/mtr_lf.csv -append 0;
sct_extract_metric -i mtr.nii.gz -f label/atlas -method map -l 55 -vert 2:4 -vertfile label/template/PAM50_levels.nii.gz -z 3:18 -o csv/mtr_vf.csv -append 0;

# Display end time

echo "${SUBJECT}: spine.sh finished at $(date +%x_%r)";
end=$(date +%s);
echo "${SUBJECT}: spine.sh duration: $((end-start)) seconds"
