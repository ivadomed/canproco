"""
This script crops the images and their label (if it exists) to the same size as the mask given. 
The mask given is the segmentation of the spinal cord.

Args:
    input_folder: path to the folder containing the images
    contrast: contrast used for the segmentation
    qc_folder: path to the quality control folder
    replacing: if True, the images will be replaced by the cropped images and the segmentations will be removed. If False, the cropped images will be saved with the suffix '_crop'.

Returns:
    None

Example:
    python crop_dataset.py -i /path/to/dataset -c PSIR,STIR -q /path/to/qc/folder

Todo:
    *

Pierre-Louis Benveniste
"""

import os
import argparse
from pathlib import Path

def crop_to_mask(image, mask, output_path, qc_folder, dilate_size='2x1000x2',label=False):
    """
    This function is used to crop the images at the same size as the mask given. 
    It dilates the mask to avoid cutting the edges of the mask.
    It also creates the quality control on the sagittal plane.

    Args:
        image: image to crop
        mask: mask to crop
        dilate_size: size of the dilation
        output_path: output folder path
    
    Returns:
        None
    """
    
    os.system(f'sct_crop_image -i {image} -m {mask} -o {output_path} -dilate {dilate_size} -b 0')

    if not label:
        os.system(f'sct_qc -i {image} -d {output_path} -p sct_image_stitch -plane sagittal -qc {qc_folder}')

    return None


def replacing_fc(image, mask, cropped,label=False):
    """
    This function replaces the image by the cropped image and removes the segmentation.
    It also deleted the _centerline.nii.gz file if it exists. 


    Args:
        image: image to replace
        mask: mask to replace
        cropped: cropped image
    
    Returns:
        None
    """
    if label:
        os.system(f'rm {image}')
        os.system(f'mv {cropped} {image}')

    else:
        os.system(f'rm {image}')
        os.system(f'rm {mask}')
        os.system(f'mv {cropped} {image}')

        if os.path.exists(str(image).split('.')[0] + '_centerline.nii.gz'):
            os.system(f'rm {str(image).split(".")[0] + "_centerline.nii.gz"}')

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
    parser.add_argument('-q', '--qc_folder', type=str, help='path to the quality control folder', required=True)
    parser.add_argument('-r', '--replacing', type=bool, help='if True, the images will be replaced by the cropped images and the segmentations will be removed. If False, the cropped images will be saved with the suffix _crop', required=False, default=False)
    
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
        crop_to_mask(image, mask, output_image, args.qc_folder)

        #get the label
        label = str(input_folder) + '/derivatives/labels/' + mouse_nb  + '/' + session_nb + '/anat/' + str(image).split('/')[-1].split('.')[0] + '_lesion-manual.nii.gz'
        #if the label exists
        if os.path.exists(label):
            #output label
            output_label = label.split('.')[0] + '_crop.nii.gz'
            #Crop the label
            crop_to_mask(label, mask, output_label, args.qc_folder,label=True)
            #if replacing is True
            if args.replacing:
                #replace the label by the cropped label
                replacing_fc(label, mask, output_label,label=True)
        
        #if replacing is True
        if args.replacing:
            #replace the image by the cropped image
            replacing_fc(image, mask, output_image)

    return None

if __name__ == '__main__':
    main()