"""
This function is used to crop the images at the same size as the mask given

Args:
    image: image to crop
    mask: mask to crop
    dilate_size: size of the dilation
    output_path: output folder path

Returns:
    None

Todo:
    * 

Pierre-Louis Benveniste
"""

import os

def crop_to_mask(image, mask, output_path, dilate_size=2):
    #This function uses sct_crop_image to crop the image to the mask
    #The mask is dilated to avoid cutting the edges of the mask
    #The image is cropped to the mask

    os.system(f'sct_crop_image -i {image} -m {mask} -o {output_path} -dilate {dilate_size} ')

    return None
