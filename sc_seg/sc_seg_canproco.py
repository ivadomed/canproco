"""
This scripts segments the spinal cord of every subject in the CANPROCO dataset using the contrast agnostic model. 
For PSIR images, the image is multiplied by -1 before segmentation as it improves the quality of the segmentation.

Args:
    input_data : path to the BIDS dataset
    contrast: contrasts to crop (separated by commas)
    qc_folder: path to the quality control folder
    seg_script: path to the script that segments the spinal cord
    path_to_model: path to the model that segments the spinal cord
    output_folder: path to the output folder
    time_point: time point to segment (optional)

Returns:
    None

Example:
    python sc_seg_canproco.py --input_data /path/to/dataset --contrast PSIR,STIR --qc_folder /path/to/qcFoldert --seg_script /path/to/script --path_to_model /path/to/model --output_folder /path/to/output/folder

Todo:
    *

Pierre-Louis Benveniste
"""


import os
import argparse
from pathlib import Path


def get_parser():
    """
    This function parses the command line arguments and returns an argparse object.

    Input:
        None

    Returns:
        parser : argparse object
    """
    parser = argparse.ArgumentParser(description='Segments the spinal cord of every subject in the CANPROCO dataset using the contrast agnostic model.')
    parser.add_argument('--input-data', '-i', type=str, help='path to the BIDS dataset')
    parser.add_argument('--contrast', '-c', type=str, help='contrasts to crop (separated by commas)')
    parser.add_argument('--qc-folder', '-q', type=str, help='path to the quality control folder')
    parser.add_argument('--seg-script', '-s', type=str, help='path to the script that segments the spinal cord')
    parser.add_argument('--path-to-model', '-m', type=str, help='path to the model that segments the spinal cord')
    parser.add_argument('--output-folder', '-o', type=str, help='path to the output folder')
    parser.add_argument('--time-point', '-t', type=str, help='time point to segment (optional)', default='M0')

    return parser


def binarize_mask(mask_path):
    """This function takes a mask and binarizes it using sct_maths.

    Args:
        mask_path: path to the mask

    Returns:
        None
    """

    os.system(f'sct_maths -i {mask_path} -o {mask_path} -bin 0.5 ')

    return None



def sc_seg_subject(img_path, qc_folder, path_to_seg_script, path_to_model, output_folder):
    """
    This function performs the segmentation of the spinal cord of a subject using the contrast agnostic model.
    If the contrast is PSIR, the image is multiplied by -1 before segmentation.

    Args:
        img_path: path to the subject
        qc_folder: path to the quality control folder
        path_to_seg_script: path to the script that segments the spinal cord
        path_to_model: path to the model that segments the spinal cord
        output_folder: path to the output folder
    
    Returns:
        output_file : path to the segmentation of the spinal cord
    """
    #create temporary folder
    tmp_folder = os.path.join(output_folder, 'tmp')
    if not os.path.exists(tmp_folder):
        os.mkdir(tmp_folder)

    

    #multiply the image by -1 if contrast is PSIR
    if 'PSIR' in str(img_path):
        #create temp file
        tmp_file = os.path.join(tmp_folder, str(img_path).split('/')[-1])
        os.system(f'sct_maths -i {img_path} -o {tmp_file} -mul -1 ')
        #segment the spinal cord
        os.system(f'python {path_to_seg_script} --path-img {tmp_file}  --chkp-path {path_to_model} --path-out {output_folder} --crop-size 100x320x320 --device gpu ')

    else:
        #segment the spinal cord
        os.system(f'python {path_to_seg_script} --path-img {img_path}  --chkp-path {path_to_model} --path-out {output_folder} --crop-size 100x320x320 --device gpu ')

    output_file = os.path.join(output_folder,str(img_path).split('/')[-1].replace('.nii.gz', '_pred.nii.gz'))
    
    #binarize the mask
    binarize_mask(output_file)

    #create qc
    os.system(f'sct_qc -i {img_path} -s {output_file} -d {output_file} -p sct_deepseg_lesion -plane sagittal -qc {qc_folder}')

    return output_file


def main():
    """
    This function performs the segmentation of the spinal cord of every subject in the CANPROCO dataset using the contrast agnostic model.

    Args: 
        None

    Returns:
        None
    """

    #Get the arguments
    parser = get_parser()
    args = parser.parse_args()
    input_data = Path(args.input_data)
    contrasts = list(args.contrast.split(','))
    qc_folder = Path(args.qc_folder)
    seg_script = Path(args.seg_script)
    path_to_model = Path(args.path_to_model)
    output_folder = Path(args.output_folder)

    #Create the output folder if it does not exist
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    
    #Loop over the subjects with desired contrasts
    image_files = []
    for contrast in contrasts:
        image_files += list(input_data.rglob(f'*_{contrast}.nii.gz'))
    image_files = sorted(image_files)
    
    #keep those only of wanted timepoint
    image_files = [image_file for image_file in image_files if args.time_point in str(image_file)]

    #loop over the image file
    for image_file in image_files:
        #segment the spinal cord using the contrast agnostic model
        output_file = sc_seg_subject(image_file, qc_folder, seg_script, path_to_model, output_folder)


    return None


if __name__ == '__main__':
    main()