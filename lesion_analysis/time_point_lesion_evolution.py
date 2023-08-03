"""
This python file is used to analyse the evolution of the lesions over time.
It first performs a registration of the images and their segmentation on the two time points.
Then, it performs clustering of the lesions on each image and computes the volume and the center of each lesion. 
Finally, it compares the lesions between the two time points and computes the evolution of the lesions.

Args:
    -i-1, --input-first_image: path to the first image 
    -i-2, --input-second_image: path to the second image
    -seg-1, --segmentation-first_image: path to the lesion segmentation of the first image
    -seg-2, --segmentation-second_image: path to the lesion segmentation of the second image
    -o, --output_folder: path to the output folder
    --plot: whether to plot the results or not

Returns:
    None

Example:
    python time_point_lesion_evolution.py -i-1 path/to/first/image -i-2 path/to/second/image -seg-1 path/to/first/segmentation -seg-2 path/to/second/segmentation -o path/to/output/folder
    python3 lesion_analysis/time_point_lesion_evolution.py -i-1 /Users/plbenveniste/Desktop/lesion_comparison_copy/sub-cal072_ses-M0_STIR.nii.gz -i-2 /Users/plbenveniste/Desktop/lesion_comparison_copy/sub-cal072_ses-M12_STIR.nii.gz -seg-1 /Users/plbenveniste/Desktop/lesion_comparison_copy/M0_inference_results.nii.gz -seg-2 /Users/plbenveniste/Desktop/lesion_comparison_copy/M12_inference_results.nii.gz -o /Users/plbenveniste/Desktop/lesion_comparison_copy/output

To do:
    *

Pierre-Louis Benveniste
"""

import os
import argparse
from pathlib import Path
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN

def get_parser():
    """
    This function parses the arguments given to the script.

    Args:
        None

    Returns:
        parser: parser containing the arguments
    """

    parser = argparse.ArgumentParser(description='Analyse the results of the segmentation of the MS lesion on the spinal cord.')
    parser.add_argument('-i-1', '--input-first_image', required=True,
                        help='Path to the first image')
    parser.add_argument('-i-2', '--input-second_image', required=True,
                        help='Path to the second image')
    parser.add_argument('-seg-1', '--segmentation-first_image', required=True,
                        help='Path to the lesion segmentation of the first image')
    parser.add_argument('-seg-2', '--segmentation-second_image', required=True,
                        help='Path to the lesion segmentation of the second image')
    parser.add_argument('-o', '--output_folder', required=True,
                        help='Path to the output folder')
    parser.add_argument('--show', action='store_true',
                        help='Whether to show the results or not')
    return parser

def main():
    """
    This function is the main function of the script.
    
    Args:
        None

    Returns:
        None
    """
    #get the parser
    parser = get_parser()
    args = parser.parse_args()

    #get the images and the segmentations
    first_image = nib.load(args.input_first_image)
    second_image = nib.load(args.input_second_image)
    first_lesion_segmentation = nib.load(args.segmentation_first_image)
    second_lesion_segmentation = nib.load(args.segmentation_second_image)

    #first we segment the spinal cord on the two images using sct_deepseg_sc
    #os.system('sct_deepseg_sc -i ' + args.input_first_image + ' -c t2 -o ' + args.output_folder + '/first_image_sc_segmentation.nii.gz')
    #os.system('sct_deepseg_sc -i ' + args.input_second_image + ' -c t2 -o' + args.output_folder + '/second_image_sc_segmentation.nii.gz')

    #then we compute the vertebral levels of the images using sct_label_vertebrae
    #os.system('sct_label_vertebrae -i ' + args.input_first_image + ' -s ' + args.output_folder + '/first_image_sc_segmentation.nii.gz' + ' -c t2 -qc ' + args.output_folder + '/qc' + ' -ofolder ' + args.output_folder)
    #os.system('sct_label_vertebrae -i ' + args.input_second_image + ' -s ' + args.output_folder + '/second_image_sc_segmentation.nii.gz' + ' -c t2 -qc ' + args.output_folder + '/qc' + ' -ofolder ' + args.output_folder)

    
    #then we register the images and the segmentations with the first image as reference using sct_register_multimodal using both the spinal cord segmentation and the vertebral levels
    parameters = 'step=0,type=label,dof=Tx_Ty_Tz:step=1,type=seg,algo=slicereg,poly=5:step=2,type=im,algo=syn,metric=MI,deformation=1x1x1,smooth=3'
    os.system('sct_register_multimodal -i ' + args.input_second_image + ' -d ' + args.input_first_image + ' -iseg ' + args.output_folder + '/second_image_sc_segmentation.nii.gz'
                + ' -dseg ' + args.output_folder + '/first_image_sc_segmentation.nii.gz' + ' -ilabel ' + args.output_folder + '/second_image_sc_segmentation_labeled.nii.gz' 
                + ' -dlabel ' + args.output_folder + '/first_image_sc_segmentation_labeled.nii.gz' + ' -o ' + args.output_folder + '/first_image_registered.nii.gz'
                + ' -owarp ' + args.output_folder +'/warping_M0_to_M0.nii.gz' + ' -param ' + parameters + ' -x linear -qc ' + args.output_folder + '/qc')
                

    # parameters = 'step=1,type=seg,algo=slicereg,poly=5:step=2,type=im,algo=syn,metric=MI,deformation=1x1x1,smooth=3'
    # os.system('sct_register_multimodal -i ' + args.input_second_image + ' -d ' + args.input_first_image + ' -iseg ' + args.output_folder + '/second_image_sc_segmentation.nii.gz' 
    #            + ' -dseg ' + args.output_folder + '/first_image_sc_segmentation.nii.gz' + ' -o ' + args.output_folder + '/second_image_registered.nii.gz' 
    #            + ' -owarp ' + args.output_folder +'/warping_M0_to_M12.nii.gz' + ' -param ' + parameters + ' -x linear -qc ' + args.output_folder + '/qc')
    # #then we apply the warping field to the segmentation of the lesions on the second images using sct_apply_transfo
    # os.system('sct_apply_transfo -i ' + args.segmentation_second_image + ' -d ' + args.segmentation_first_image + ' -w ' + args.output_folder + '/warping_M0_to_M12.nii.gz' + ' -o ' + args.output_folder + '/second_image_lesion_segmentation_registered.nii.gz' + ' -x linear')

    return None

if __name__ == '__main__':
    main()