"""
This python file is used to analyse the evolution of the lesions over time.
It first performs a registration of the images and their segmentation on the two time points.
Then, it performs clustering of the lesions on each image and computes the volume and the center of each lesion. 
Finally, it compares the lesions between the two time points and computes the evolution of the lesions.

Args:
    -i-1, --input-first_image: path to the first image 
    -i-2, --input-second_image: path to the second image
    -seg-1, --segmentation-first_image: path to the segmentation of the first image
    -seg-2, --segmentation-second_image: path to the segmentation of the second image
    -o, --output_folder: path to the output folder
    --plot: whether to plot the results or not

Returns:
    None

Example:
    python time_point_lesion_evolution.py -i-1 path/to/first/image -i-2 path/to/second/image -seg-1 path/to/first/segmentation -seg-2 path/to/second/segmentation -o path/to/output/folder

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
                        help='Path to the segmentation of the first image')
    parser.add_argument('-seg-2', '--segmentation-second_image', required=True,
                        help='Path to the segmentation of the second image')
    parser.add_argument('-o', '--output_folder', required=True,
                        help='Path to the output folder')
    parser.add_argument('--plot', action='store_true',
                        help='Whether to plot the results or not')
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
    first_segmentation = nib.load(args.segmentation_first_image)
    second_segmentation = nib.load(args.segmentation_second_image)

    #first we perform a registration of the images and the segmentations
    #we use the first image and the first segmentation as the reference
    

    return None

if __name__ == '__main__':
    main()