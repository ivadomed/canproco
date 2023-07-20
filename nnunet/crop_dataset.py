"""
This script crops the images and their label (if it exists) to the same size as the mask given. 
The mask given is the segmentation of the spinal cord.

Args:
    input_folder: path to the folder containing the images
    contrast: contrast used for the segmentation

Returns:
    None

Todo:
    *

Pierre-Louis Benveniste
"""

import os
import argparse
from pathlib import Path

def crop_to_mask(image, mask, output_path, dilate_size=2):
    """
    This function is used to crop the images at the same size as the mask given. 
    It dilates the mask to avoid cutting the edges of the mask.

    Args:
        image: image to crop
        mask: mask to crop
        dilate_size: size of the dilation
        output_path: output folder path
    
    Returns:
        None
    """
    
    os.system(f'sct_crop_image -i {image} -m {mask} -o {output_path} -dilate {dilate_size} ')
    
    return None


def get_parser():
    """
    This function parses the arguments given to the script

    Args:
        None
    
    Returns:
        parser: parser containing the arguments
    """

    parser = argparse.ArgumentParser(description='This script crops the images and their label to the same size as the mask given. The mask given is the segmentation of the spinal cord.')
    parser.add_argument('-i', '--input_folder', type=str, help='path to the folder containing the images', required=True)
    parser.add_argument('-c', '--contrast', type=str, help='contrast used for the segmentation (if multiple: do the following: PSIR,STIR (no space))', required=True)
    
    return parser

def main():
    """
    This function is the main function of the script. It calls the other functions to crop the images and its label (if it exists).
    It performs the cropping for all the images in the BIDS format dataset

    Args:
        None
    
    Returns:
        None
    """

    #Get the arguments
    parser = get_parser()
    args = parser.parse_args()
    input_folder = Path(args.input_folder)
    contrasts = list(args.contrast.split(','))

    #Get the list of images
    image_files = []
    for contrast in contrasts:
        image_files += list(input_folder.rglob(f'*{contrast}.nii.gz'))

    #Crop the images and their label
    for image in image_files:
        #get the mouse number
        mouse_nb = str(image.parts[-4])
        #get session number
        session_nb = str(image.parts[-3])
        #Get the mask
        mask = str(image).split('.')[0] + '_seg.nii.gz'
        #output image
        output_image = str(image).split('.')[0] + '_crop.nii.gz'
        #Crop the image
        crop_to_mask(image, mask, output_image)
        
        #get the label
        label = str(input_folder) + '/derivatives/labels/' + mouse_nb  + '/' + session_nb + '/anat/' + str(image).split('/')[-1].split('.')[0] + '_lesion-manual.nii.gz'
        #if the label exists
        if os.path.exists(label):
            #output label
            output_label = label.split('.')[0] + '_crop.nii.gz'
            #Crop the label
            crop_to_mask(label, mask, output_label)


    return None

if __name__ == '__main__':
    main()