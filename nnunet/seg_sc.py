"""
This functions is used to segment the spinal cord using the spinal cord toolbox

Args:
    image: image to segment
    output_path: output folder path
    contrast: contrast used for the segmentation

Returns:
    None

Todo:
    * 

Pierre-Louis Benveniste
"""

import os

def segment_sc(image, output_path, contrast='t2'):
    #This function uses sct_deepseg_sc to segment the spinal cord
    #The image is segmented and the segmentation is saved in the output folder

    os.system(f'sct_deepseg_sc -i {image} -o {output_path} -c {contrast}')

    return None
