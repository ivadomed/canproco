"""
Convert region-based nnUNet predictions to BIDS format
This scripts takes the predictions of nnUNet and converts them to BIDS format with the correct file names. 

Example:
    python convert_predictions_to_BIDS.py --pred-folder /path/to/predictions --out-folder /path/to/output/folder --conversion-dict /path/to/conversion/dict

Args:
    --pred-folder : path to the folder containing the predictions of nnUNet
    --out-folder : path to the folder where the converted predictions should be saved
    --conversion-dict : path to the conversion dictionary
    --not-imageTs : if the predictions are not in the imagesTs folder, set this flag to True

Returns:
    None

TODO:
    *

Pierre-Louis Benveniste
"""

import os
import argparse
import shutil
import json
import numpy as np
import pathlib


def get_parser():
    """
    This function parses the command line arguments and returns an argparse object.

    Input:
        None

    Returns:
        parser : argparse object
    """
    parser = argparse.ArgumentParser(description="Convert region-based nnUNet predictions to BIDS format")
    parser.add_argument("--pred-folder", type=str, required=True, help="path to the folder containing the predictions of nnUNet")
    parser.add_argument("--out-folder", type=str, required=True, help="path to the folder where the converted predictions should be saved")
    parser.add_argument("--conversion-dict", type=str, required=True, help="path to the conversion dictionary")
    parser.add_argument("--not-imageTs", action="store_true", help="if the predictions are not in the imagesTs folder, set this flag to True")
    return parser


def main():
    """
    This function executes the conversion of the nnUNet predictions to BIDS format.

    Input:
        None

    Returns:
        None
    """

    #we parse the command line arguments
    parser = get_parser()
    args = parser.parse_args()
    pred_folder = args.pred_folder
    out_folder = args.out_folder
    conversion_dict = args.conversion_dict

    #we build the output folder if it does not exist
    if not os.path.isdir(out_folder):
        os.makedirs(out_folder)
    
    #we load the conversion dictionary
    with open(conversion_dict, "r") as f:
        conversion_dict = json.load(f)

    # we list all the files in the predictions which end with .nii.gz
    pred_files = list(pathlib.Path(pred_folder).rglob(f'*.nii.gz'))

    # we loop over the files
    for pred in pred_files:

        # we get the file name
        pred_name = pred.name
        
        # we get the corresponding file in the conversion dictionary
        for files in conversion_dict:

            if (pred_name.split('.')[0] in conversion_dict[files].split("/")[-1] and 'imagesTs' in conversion_dict[files] ) or (pred_name.split('.')[0] in conversion_dict[files].split("/")[-1] and args.not_imageTs) :
                
                subject_full = files.split("/")[-1]
                subject = subject_full.split("_")[0]
                session = subject_full.split("_")[1]

                #we create the corresponding folder
                output_subject_path = os.path.join(out_folder, subject, session, 'anat') 
                if not os.path.isdir(output_subject_path):
                    os.makedirs(output_subject_path)
                # and the new name
                new_name = subject_full.split('.')[0] + '_pred.nii.gz'
                
                #we copy the file to the corresponding folder and rename it
                shutil.copy(pred, os.path.join(output_subject_path, new_name))
    
    return None


if __name__ == "__main__":
    main()