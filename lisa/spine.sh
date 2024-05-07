#!/bin/bash
# Written by Lisa Eunyoung Lee - Jun 2022
# Revised by Lisa Eunyoung Lee - Jan 2023
# Spinal Cord Toolbox - Version 5.7

# T2 - Spinal Cord CSA averaged between C2-C4
# MT - Spinal Cord MTR averaged between C2-C4

# Check input arguments

if [ $# -eq 0 ] ; then
    echo "spine.sh"
	echo "sh spine.sh SITE GROUP NUMSUBS MONTH"
	printf "Options: \n SITE = 01 (SMH Toronto), 02 (UBC Vancouver), 03 (CHUM Montreal), 04 (UCAL Calgary), 05 (UALB Edmonton) \n GROUP = RIS (RIS), RRM (RRMS), PPM (PPMS), CON (control) \n NUMSUBS = subject number (e.g. '2 56 105') \n MONTH = 0 (baseline), 12 (month 12), 24 (month 24), 48 (month 48) \n"
	exit 1
elif [ $# -ne 4 ] ; then
	echo "spine.sh: Incorrect number of inputs, $# provided but 4 required." 1>&2
	exit 1
fi; 

# Set up variables

SITE=$1
if [ $SITE = '01' ]; then STUDYSITE='Toronto'
elif [ $SITE = '02' ]; then STUDYSITE='Vancouver'
elif [ $SITE = '03' ]; then STUDYSITE='Montreal'
elif [ $SITE = '04' ]; then STUDYSITE='Calgary'
elif [ $SITE = '05' ]; then STUDYSITE='Edmonton'
else echo "Site type not recognized. Only '01', '02', '03', '04' or '05' are accepted." 1>&2;exit 1
fi
GROUP=$2
NUMSUBS=$3
MONTH=$4

for NUMSUB in $NUMSUBS;do
	if [ $SITE = '01' ] || [ $SITE = '02' ] || [ $SITE = '03' ] || [ $SITE = '04' ] || [ $SITE = '05' ]; then NUM=$(printf "%03d\n" $NUMSUB)
	fi

# Set up paths

	data="/Volumes/Neuron/Lisa/canproco/data/spine/$STUDYSITE"
	ubc="/Volumes/Neuron/Lisa/canproco/data/ubc/$STUDYSITE/CAN-$SITE-$GROUP-$NUM-M${MONTH}/cache"
	qc="/Volumes/Neuron/Lisa/canproco/data/spine/qc_SCT/M${MONTH}/$STUDYSITE"

# Display start time

	echo "CAN-$SITE-$GROUP-$NUM-M${MONTH}: spine.sh started at $(date +%x_%r)";
	start=$(date +%s);


# Organize files
# ===========================================================================================

# 1. Organize files in the working directory

	mkdir -p "$data/CAN-$SITE-$GROUP-$NUM-M${MONTH}";
	mkdir -p "$data/CAN-$SITE-$GROUP-$NUM-M${MONTH}/t1";
	mkdir -p "$data/CAN-$SITE-$GROUP-$NUM-M${MONTH}/t1/raw";
	mkdir -p "$data/CAN-$SITE-$GROUP-$NUM-M${MONTH}/t2";
	mkdir -p "$data/CAN-$SITE-$GROUP-$NUM-M${MONTH}/t2/raw";
	mkdir -p "$data/CAN-$SITE-$GROUP-$NUM-M${MONTH}/t2s";
	mkdir -p "$data/CAN-$SITE-$GROUP-$NUM-M${MONTH}/t2s/raw";
	mkdir -p "$data/CAN-$SITE-$GROUP-$NUM-M${MONTH}/mt";
	mkdir -p "$data/CAN-$SITE-$GROUP-$NUM-M${MONTH}/mt/raw";

	data="/Volumes/Neuron/Lisa/canproco/data/spine/$STUDYSITE/CAN-$SITE-$GROUP-$NUM-M${MONTH}"

	cp $ubc/*T2W*-S.nii.gz $data/t2/raw;
	cp $ubc/*Multi_Echo*-S.nii.gz $data/t2s/raw;
	cp $ubc/*PSIR*-S.nii.gz $data/t1/raw;
	cp $ubc/json/*pine*MT*.json $data/mt/raw;
	cp $ubc/json/*pine*T1*.json $data/mt/raw;

	if [ $SITE = '02' ] || [ $SITE = '03' ]; then
		if [ -f $ubc/*MT*-S_*cale*.nii.gz -a -f $ubc/*T1*-S_*cale*.nii.gz ]; then 
			cp $ubc/*MT*-S_*cale*.nii.gz $data/mt/raw;
		else cp $ubc/*MT*-S.nii.gz $data/mt/raw;
		fi
	fi

	if [ $SITE = '01' ] || [ $SITE = '05' ]; then
		cp $ubc/*MT*-S*.nii.gz $data/mt/raw;
	fi

	sct_convert -i $data/t2/raw/*T2W*-S.nii.gz -o $data/t2/t2.nii.gz;
	sct_convert -i $data/t2s/raw/*Multi_Echo*-S.nii.gz -o $data/t2s/t2s.nii.gz;

	if [ $SITE = '02' ] || [ $SITE = '03' ]; then
		cd $data/mt/raw/;
		if [ -f *ON_OFF*cale*.nii.gz ]; then
			fslsplit *ON_OFF*cale*.nii.gz mt -t;
		else fslsplit *ON_OFF-S.nii.gz mt -t;
		fi
		mv mt0000.nii.gz mt0.nii.gz;
		mv mt0001.nii.gz mt1.nii.gz;
		mv mt*.nii.gz ../; 
		if [ -f *T1*cale*.nii.gz ]; then
			sct_convert -i *T1*cale*.nii.gz -o ../mt_t1.nii.gz;
		else sct_convert -i *T1-S.nii.gz -o ../mt_t1.nii.gz;
		fi
		cp *MT*.json ../mt0_reg.json;
		cp *MT*.json ../mt1_reg.json;
		cp *T1*.json ../mt_t1.json;
		sct_convert -i $data/t1/raw/*real*.nii.gz -o $data/t1/t1.nii.gz;
	fi

	if [ $SITE = '01' ] || [ $SITE = '05' ]; then
		cd $data/mt/raw/;
		sct_convert -i *OFF*.nii.gz -o ../mt0.nii.gz;
		sct_convert -i *ON*.nii.gz -o ../mt1.nii.gz;
		sct_convert -i *T1*.nii.gz -o ../mt_t1.nii.gz;
		cp *on*.json ../mt1_reg.json;
		cp *off*.json ../mt0_reg.json;
		cp *T1*.json ../mt_t1.json;
		sct_convert -i $data/t1/raw/PSIR-S.nii.gz -o $data/t1/t1.nii.gz;
	fi


# T2
# ===========================================================================================

# 2. Segment spinal cord

	cd $data/t2;

	sct_image -i t2.nii.gz -setorient RPI -o t2.nii.gz;

	sct_deepseg_sc -i t2.nii.gz -c t2 -qc $qc;

# 3. Label vertebral levels

	sct_label_vertebrae -i t2.nii.gz -s t2_seg.nii.gz -c t2 -qc $qc;

	sct_label_utils -i t2_seg_labeled.nii.gz -vert-body 3,5 -o t2_labels_vert.nii.gz;

	sct_qc -i t2.nii.gz -s t2_labels_vert.nii.gz -p sct_label_utils -qc $qc;

# 4. Register t2 to template space then, warp template to t2	

	sct_register_to_template -i t2.nii.gz -s t2_seg.nii.gz -l t2_labels_vert.nii.gz -c t2 -qc $qc;

	sct_warp_template -d t2.nii.gz -w warp_template2anat.nii.gz -a 0 -qc $qc;

# 5. Compute cross-sectional area averaged between C2-C4

	sct_process_segmentation -i t2_seg.nii.gz -vert 2:4 -vertfile ./label/template/PAM50_levels.nii.gz -o csa_c2c4.csv -append 0;


# MT
# ===========================================================================================

# 6. Segment spinal cord

	cd $data/mt;

	sct_image -i mt_t1.nii.gz -setorient RPI -o mt_t1.nii.gz;

	sct_deepseg_sc -i mt_t1.nii.gz -c t1 -qc $qc;

# 7. Create a mask

	sct_create_mask -i mt_t1.nii.gz -p centerline,mt_t1_seg.nii.gz -size 35mm -f cylinder -o mt_t1_mask.nii.gz;

# 8. Register template -> mt_t1 then, warp template to mt_t1 

	sct_register_multimodal -i $SCT_DIR/data/PAM50/template/PAM50_t1.nii.gz -iseg $SCT_DIR/data/PAM50/template/PAM50_cord.nii.gz -d mt_t1.nii.gz -dseg mt_t1_seg.nii.gz -m mt_t1_mask.nii.gz -param step=1,type=seg,algo=slicereg,metric=MeanSquares,smooth=2:step=2,type=im,algo=syn,metric=CC,iter=5,gradStep=0.5 -initwarp $data/t2/warp_template2anat.nii.gz -initwarpinv $data/t2/warp_anat2template.nii.gz -owarp warp_template2mt.nii.gz -qc $qc;

	sct_warp_template -d mt_t1.nii.gz -w warp_template2mt.nii.gz -a 1 -qc $qc;

# 9. Register mt0 -> mt_t1

	sct_register_multimodal -i mt0.nii.gz -d mt_t1.nii.gz -dseg mt_t1_seg.nii.gz -m mt_t1_mask.nii.gz -param step=1,type=im,algo=slicereg,metric=CC -x spline -o mt0_reg.nii.gz -qc $qc;

# 10. Register mt1 -> mt_t1

	sct_register_multimodal -i mt1.nii.gz -d mt_t1.nii.gz -dseg mt_t1_seg.nii.gz -m mt_t1_mask.nii.gz -param step=1,type=im,algo=slicereg,metric=CC -x spline -o mt1_reg.nii.gz -qc $qc;

# 11. Compute MTR

	sct_compute_mtr -mt0 mt0_reg.nii.gz -mt1 mt1_reg.nii.gz;

# 12. Extract MTR in whole cord, GM, WM, WM regions (DC, LF, VF, CST) averaged between C2-C4

	mkdir -p $data/mt/csv;

	sct_extract_metric -i mtr.nii.gz -f label/atlas -method map -l 50 -vert 2:4 -vertfile label/template/PAM50_levels.nii.gz -z 3:18 -o csv/mtr_cord.csv -append 0;
	sct_extract_metric -i mtr.nii.gz -f label/atlas -method map -l 52 -vert 2:4 -vertfile label/template/PAM50_levels.nii.gz -z 3:18 -o csv/mtr_gm.csv -append 0;
	sct_extract_metric -i mtr.nii.gz -f label/atlas -method map -l 51 -vert 2:4 -vertfile label/template/PAM50_levels.nii.gz -z 3:18 -o csv/mtr_wm.csv -append 0;
	sct_extract_metric -i mtr.nii.gz -f label/atlas -method map -l 53 -vert 2:4 -vertfile label/template/PAM50_levels.nii.gz -z 3:18 -o csv/mtr_dc.csv -append 0;
	sct_extract_metric -i mtr.nii.gz -f label/atlas -method map -l 54 -vert 2:4 -vertfile label/template/PAM50_levels.nii.gz -z 3:18 -o csv/mtr_lf.csv -append 0;
	sct_extract_metric -i mtr.nii.gz -f label/atlas -method map -l 55 -vert 2:4 -vertfile label/template/PAM50_levels.nii.gz -z 3:18 -o csv/mtr_vf.csv -append 0;
	sct_extract_metric -i mtr.nii.gz -f label/atlas -method map -l 56 -vert 2:4 -vertfile label/template/PAM50_levels.nii.gz -z 3:18 -o csv/mtr_cst.csv -append 0;

# Display end time

	echo "CAN-$SITE-$GROUP-$NUM-M${MONTH}: spine.sh finished at $(date +%x_%r)";
	end=$(date +%s);
	echo "CAN-$SITE-$GROUP-$NUM-M${MONTH}: spine.sh duration: $((end-start)) seconds";
	
done