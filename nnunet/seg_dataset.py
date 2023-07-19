"""
This script adds segmentation of the spinal cord to the dataset using the spinal cord toolbox. 

Args:
    input_folder: path to the folder containing the images
    contrast: contrast used for the segmentation
    qc_folder: path to the quality control folder

Returns:
    None

Todo:
    *

Pierre-Louis Benveniste
"""

import os
import argparse
from pathlib import Path

def seg_sc(image, contrast, output_path, qc_folder):
    """
    This function is used to segment the spinal cord using the spinal cord toolbox. Is

    Args:
        image: image to segment
        contrast: contrast used for the segmentation
        output_path: output folder path

    Returns:
        None
    """

    contrast_dict = {'STIR': 't2', 'PSIR': 't1', 'T2star':'t2s', 'T2w': 't2', 'T1w': 't1', 'MT-T1': 't1', 'MTon': 't2s' }

    if contrast=='PSIR':
        os.system(f'sct_propseg_sc -i {image} -o {output_path} -c {contrast_dict[contrast]} -qc {qc_folder}')
    else:
        os.system(f'sct_deepseg_sc -i {image} -o {output_path} -c {contrast_dict[contrast]} -qc {qc_folder}')

    return None

def get_parser():
    """
    This function parses the arguments given to the script

    Args:
        None
    
    Returns:
        parser: parser containing the arguments
    """

    parser = argparse.ArgumentParser(description='This script adds segmentation of the spinal cord to the dataset using the spinal cord toolbox.')
    parser.add_argument('-i', '--input_folder', type=str, help='path to the folder containing the images', required=True)
    parser.add_argument('-c', '--contrast', type=str, help='contrast used for the segmentation (if multiple: do the following: PSIR,STIR (no space))', required=True)
    parser.add_argument('-q', '--qc_folder', type=str, help='path to the quality control folder', required=True)
    
    return parser

def main():
    """
    This function is the main function of the script. It calls the other functions to segment the spinal cord.

    Args:
        None

    Returns:
        None
    """

    parser = get_parser()
    args = parser.parse_args()

    #get the path to the folder containing the images
    path_in_images = Path(args.input_folder)
    #get the contrast used for the segmentation
    contrasts = list(args.contrast.split(','))

    #Get the list of segmentations with the right contrast
    seg_files = []
    for contrast in contrasts:
        seg_files += sorted(list(path_in_images.rglob(f'*{contrast}.nii.gz')))
    seg_files = sorted(list(set(seg_files)))

    for seg_file in seg_files:
        #get the mouse number
        mouse_nb = seg_file.parts[-4]
        #get session number
        session_nb = seg_file.parts[-3]
        #get the contrast
        contrast = seg_file.parts[-1].split('.')[0].split('_')[-1]
        #get the output path
        seg_out_path = str(seg_file).split('.')[0] + '_seg.nii.gz'
        #segment the spinal cord
        seg_sc(seg_file, contrast, seg_out_path, args.qc_folder)

if __name__ == '__main__':
    main()