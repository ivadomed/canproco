"""
This script is used to evaluate the lesion segmentation prediction of the test files.
It looks at different metrics (dice score, precision and recall) and saves the results in a csv file.
We use anima for the computation of the metrics. Here is the return message from animaSegPerfAnalyzer:
------------------------------------------------------------------------------------------------------------------------
3 categories are available:
    - SEGMENTATION EVALUATION:
        Dice, the mean overlap
        Jaccard, the union overlap
        Sensitivity
        Specificity
        NPV (Negative Predictive Value)
        PPV (Positive Predictive Value)
        RVE (Relative Volume Error) in percentage
    - SURFACE DISTANCE EVALUATION:
        Hausdorff distance
        Contour mean distance
        Average surface distance
    - DETECTION LESIONS EVALUATION:
        PPVL (Positive Predictive Value for Lesions)
        SensL, Lesion detection sensitivity
        F1 Score, a F1 Score between PPVL and SensL

Results are provided as follows: 
Jaccard;        Dice;   Sensitivity;    Specificity;    PPV;    NPV;    RelativeVolumeError;    HausdorffDistance;      ContourMeanDistance;    SurfaceDistance;        PPVL;   SensL;  F1_score;       NbTestedLesions;        VolTestedLesions;
------------------------------------------------------------------------------------------------------------------------

Args:
    --i, --input_img : path to the input image
    --l, --lesion_mask : path to the lesion mask
    --p, --prediction : path to the prediction
    --animaPath : path to the animaSegPerfAnalyzer script
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
import xml.etree.ElementTree as ET


def get_parser():
    """
    This function parses the command line arguments and returns an argparse object.

    Input:
        None

    Returns:
        parser : argparse object
    """
    parser = argparse.ArgumentParser(description='Evaluate lesion segmentation prediction')
    parser.add_argument('--i', required=True, type=str, help='path to the input image')
    parser.add_argument('--l',required=True, type=str, help='path to the lesion mask')
    parser.add_argument('--p', required=True, type=str, help='path to the prediction')
    parser.add_argument('--animaPath', required=True, type=str, help='path to the animaSegPerfAnalyzer script')
    parser.add_argument('--o', required=True,  type=str, help='path to the output folder')
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
    input_img = args.i
    lesion_mask = args.l
    prediction = args.p
    animaPath = args.animaPath
    output_folder = args.o

    #get name
    name = input_img.split('/')[-1].split('.')[0]

    #we build the output folder if it does not exist
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)
    if not os.path.isdir(output_folder + '/tmp/' + name):
        os.makedirs(output_folder + '/tmp/'+ name)
    if not os.path.isdir(output_folder + '/' + name):
        os.makedirs(output_folder + '/' + name)
    
    #we binarize the prediction and the lesion mask using sct_maths
    prediction_bin_path = output_folder + '/tmp/' + name + '/' + prediction.split('/')[-1].split('.')[0] + '_bin.nii.gz'
    lesion_mask_bin_path = output_folder + '/tmp/' + name + '/' + lesion_mask.split('/')[-1].split('.')[0] + '_bin.nii.gz'
    os.system('sct_maths -i ' + prediction + ' -bin 0.5 -o ' + prediction_bin_path)
    os.system('sct_maths -i ' + lesion_mask + ' -bin 0.5 -o ' + lesion_mask_bin_path)

    #we use the script from Anima : animaSegPerfAnalyzer
    print('Computing metrics with Anima...')
    output_path = output_folder + '/tmp/' + name + '/' + name + '_anima_analysis'
    print(output_path)
    os.system( animaPath + '/animaSegPerfAnalyzer -r ' + lesion_mask_bin_path + ' -i ' + prediction_bin_path + ' -o ' + output_path + '-d -s -l -X -S')
    print('Done!')

    #we initalise the dtaframe of the subject results
    subject_results = {'subject': prediction.split('/')[-1].split('.')[0]}
    
    #we read the xml file
    xml_file = output_path + '-d_global.xml'
    root_node = ET.parse(source=xml_file).getroot()
    for metric in root_node:
        metric_name, value = metric.get('name'), float(metric.text)
        subject_results[metric_name] = value
    
    #save the results in a csv file
    df = pd.DataFrame(subject_results, index=[0])
    df.to_csv(output_folder + '/' + name + '/' + name + '_results.csv', index=False)


if __name__ == "__main__":
    main()
