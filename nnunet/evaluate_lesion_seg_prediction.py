"""
This script is used to evaluate the lesion segmentation prediction.
It looks at different metrics (dice score, precision and recall) and saves the results in a csv file.
List of metrics outputted:
- Dice score
- Precision
- Recall
- Lesion wide precision 
- Lesion wide recall
- TP, FP, TN, FN

Args:
    --i, --input_img : path to the input image
    --l, --lesion_mask : path to the lesion mask
    --p, --prediction : path to the prediction
    --o, --output_folder : path to the output folder

Returns:
    None

Example:
    python evaluate_lesion_seg_prediction.py --i /path/to/image --l /path/to/lesion/mask --p /path/to/prediction --o /path/to/output/folder

TODO: 
    *

Pierre-Louis Benveniste
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd


def get_parser():
    """
    This function parses the command line arguments and returns an argparse object.

    Input:
        None

    Returns:
        parser : argparse object
    """
    parser = argparse.ArgumentParser(description='Evaluate lesion segmentation prediction')
    parser.add_argument('--i', '--input_img', required=True, type=str, help='path to the input image')
    parser.add_argument('--l', '--lesion_mask',required=True, type=str, help='path to the lesion mask')
    parser.add_argument('--p', '--prediction', required=True, type=str, help='path to the prediction')
    parser.add_argument('--o', '--output_folder', required=True,  type=str, help='path to the output folder')
    args = parser.parse_args()
    return args


def main():
    """
    This function is the main function of the script.

    Args:
        None
    
    Returns:
        None
    """
    #we parse the command line arguments
    args = get_parser()

    #we load the image, the lesion mask and the prediction
    img = np.load(args.input_img)
    lesion_mask = np.load(args.lesion_mask)
    prediction = np.load(args.prediction)

    #we build the output folder if it does not exist
    if not os.path.isdir(args.output_folder):
        os.makedirs(args.output_folder)
    
    #we first compute the dice score
    #using anima ? anima is only runnable on Linux, not on Mac...
    
     
