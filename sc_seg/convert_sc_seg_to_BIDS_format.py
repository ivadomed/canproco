"""
This file is used to convert the segmentation of the spinal cord of every subject in the CANPROCO dataset to the BIDS format.
It copies the segmentation from a source folder to the CanProCo dataset and creates a json file for every segmentation.

Args:
    seg_folder: path to the folder containing the segmentations of the spinal cord
    bids_folder: path to the CanProCo dataset

Returns:
    None

Example:
    python convert_sc_seg_to_BIDS_format.py --seg_folder /path/to/segmentations --bids_folder /path/to/bids

Todo:
    *

Pierre-Louis Benveniste
"""


import os
import argparse
from pathlib import Path
from datetime import datetime
import json
import shutil


def get_parser():
    """
    This function parses the command line arguments and returns an argparse object.

    Input:
        None

    Returns:
        parser : argparse object
    """
    parser = argparse.ArgumentParser(description='Converts the segmentation of the spinal cord of every subject in the CANPROCO dataset to the BIDS format.')
    parser.add_argument('--seg_folder', '-s', type=str, help='path to the folder containing the segmentations of the spinal cord')
    parser.add_argument('--bids_folder', '-b', type=str, help='path to the CanProCo dataset')

    return parser


def main():
    """
    This function converts the segmentation of the spinal cord of every subject in the CANPROCO dataset to the BIDS format.

    Input:
        None

    Returns:
        None
    """
    parser = get_parser()
    args = parser.parse_args()

    seg_folder = args.seg_folder
    bids_folder = args.bids_folder

    seg_files = list(seg_folder.rglob(f'*.nii.gz'))

    for file in seg_files:
        file_name = str(file).split('/')[-1]
        subject = file_name.split('_')[0]
        time_point = file_name.split('_')[1]
        contrast = file_name.split('_')[2]

        output_file_name = file_name.split('_pred')[0] + '_sc_seg.nii.gz'
        output_file_path = os.path.join(bids_folder, 'derivatives', 'labels', subject, time_point, 'anat', output_file_name)
        shutil.copy2(file, output_file_path)

        #write A JSON file for each segmentation:
        json_file_name = file_name.split('_pred')[0] + '_sc_seg.json'
        json_file_path = os.path.join(bids_folder, 'derivatives', 'labels',subject, time_point, 'anat', json_file_name)
        json_data = {
            "Author": "Pierre-Louis Benveniste",
            "Date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        with open(json_file_path, "w") as json_file:
            json.dump(json_data, json_file, indent=4)

    return None


if __name__ == '__main__':
    main()