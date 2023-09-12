"""Convert data from BIDS to nnU-Net format
This scripts was adapted from a script by the following authors : Julian McGinnis, Naga Karthik 
This python script converts data from the BIDS format to the nnU-Net format in order to be able to perform pre-processing, training and inference.
This script is specifically designed to the canproco dataset. It converts the dataset from a BIDS format to a nnU-Net format. 
This script is specific in that it takes all labeled data for training and validation. Because, very little labelled data is available, no labelled data is used for testing.

Example of run:
    $ python convert_bids_to_nnunet.py --path-data /path/to/data_extracted --path-out /path/to/nnUNet_raw --taskname TASK-NAME --tasknumber DATASET-ID

Arguments:
    --path-data : Path to BIDS structured dataset. Accepts both cross-sectional and longitudinal datasets
    --path-out : Path to output directory.
    --taskname: Specify the task name - usually the anatomy to be segmented, e.g. Hippocampus
    --tasknumber : Specify the task number, has to be greater than 100 but less than 999
    --label-folder : Path to the label folders in derivatives (default='labels')
    --contrasts : Contrasts used for the images (default='PSIR') (separated by a comma)

Returns:
    None
    
Todo:
    * 
Pierre-Louis Benveniste
"""

import argparse
import pathlib
from pathlib import Path
import json
import os
import shutil
from collections import OrderedDict
from tqdm import tqdm

import nibabel as nib
import numpy as np


def get_parser():
    """
    This function parses the command line arguments and returns an argparse object.

    Input:
        None

    Returns:
        parser : argparse object
    """
    parser = argparse.ArgumentParser(description='Convert BIDS-structured database to nnUNet format.')
    parser.add_argument('--path-data', required=True,
                        help='Path to BIDS structured dataset. Accepts both cross-sectional and longitudinal datasets')
    parser.add_argument('--path-out', help='Path to output directory.', required=True)
    parser.add_argument('--taskname', default='MSSpineLesion', type=str,
                        help='Specify the task name - usually the anatomy to be segmented, e.g. Hippocampus',)
    parser.add_argument('--tasknumber', default=501,type=int, 
                        help='Specify the task number, has to be greater than 500 but less than 999. e.g 502')
    parser.add_argument('--label-folder', help='Path to the label folders in derivatives', default='labels', type=str)
    parser.add_argument('--contrasts', help='Contrast used for the images', default='PSIR', type=str)
    return parser


def main():
    """
    This function is the main function of the script. It converts the data from BIDS to nnU-Net format.

    Input:
        None

    Returns:
        None
    """
    #------------- PARSE COMMAND LINE ARGUMENTS --------------------------
    parser = get_parser()
    args = parser.parse_args()

    path_in_images = Path(args.path_data)
    label_folder = args.label_folder
    path_in_labels = Path(os.path.join(args.path_data, 'derivatives', label_folder))
    path_out = Path(os.path.join(os.path.abspath(args.path_out), f'Dataset{args.tasknumber}_{args.taskname}'))
    contrasts = args.contrasts.split(',')

    # define paths for train and test folders 
    path_out_imagesTr = Path(os.path.join(path_out, 'imagesTr'))
    path_out_imagesTs = Path(os.path.join(path_out, 'imagesTs'))
    path_out_labelsTr = Path(os.path.join(path_out, 'labelsTr'))
    path_out_labelsTs = Path(os.path.join(path_out, 'labelsTs'))

    # we load both train and validation set into the train images as nnunet uses cross-fold-validation
    train_images, train_labels = [], []
    test_images, test_labels = [], []

    # make the directories
    pathlib.Path(path_out).mkdir(parents=True, exist_ok=True)
    pathlib.Path(path_out_imagesTr).mkdir(parents=True, exist_ok=True)
    pathlib.Path(path_out_imagesTs).mkdir(parents=True, exist_ok=True)
    pathlib.Path(path_out_labelsTr).mkdir(parents=True, exist_ok=True)
    pathlib.Path(path_out_labelsTs).mkdir(parents=True, exist_ok=True)

    conversion_dict = {}

    #------------- EXTRACTION OF THE LABELED IMAGES NAMES--------------------------
    labelled_imgs = []
    
    # We first extract all the label files' names for every wanted contrast
    label_files = []
    for contrast in contrasts:
        label_files += list(path_in_labels.rglob(f'*_{contrast}_lesion-manual.nii.gz'))
    label_files = sorted(label_files)   
    labelled_imgs += [str(k) for k in label_files]
    
    #--------------- DISPACTH OF LABELLED IMAGES AND UNLABELED IMAGES ------------------- 
    
    #Initialise the number of scans in train and in test folder
    scan_cnt_train, scan_cnt_test = 0, 0

    valid_train_imgs = []
    valid_test_imgs = []

    #The image folders
    image_files = []
    for contrast in contrasts:
        image_files += list(path_in_images.rglob(f'*_{contrast}.nii.gz'))
    imahe_files = sorted(image_files)
    
    for image_file in image_files:

        #Identify common data path
        common = os.path.relpath(os.path.dirname(image_file), args.path_data)

        file_id = str(image_file).rsplit('_',1)[0].split('/')[-1] + "_"
        
        if (file_id in str(labelled_imgs)):
            label_file = [k for k in labelled_imgs if file_id in k][0]

            scan_cnt_train+= 1

            image_file_nnunet = os.path.join(path_out_imagesTr,f'{args.taskname}_{scan_cnt_train:03d}_0000.nii.gz')
            label_file_nnunet = os.path.join(path_out_labelsTr,f'{args.taskname}_{scan_cnt_train:03d}.nii.gz')
            
            train_images.append(str(image_file_nnunet))
            train_labels.append(str(label_file_nnunet))
            
            # copy the files to new structure
            shutil.copyfile(image_file, image_file_nnunet)
            shutil.copyfile(label_file, label_file_nnunet)
           
            conversion_dict[str(os.path.abspath(image_file))] = image_file_nnunet
            conversion_dict[str(os.path.abspath(label_file))] = label_file_nnunet
        else:
                            
            # Repeat the above procedure for testing
            scan_cnt_test += 1
            # create the new convention names
            image_file_nnunet = os.path.join(path_out_imagesTs,f'{args.taskname}_{scan_cnt_test:03d}_0000.nii.gz')
            label_file_nnunet = os.path.join(path_out_labelsTs,f'{args.taskname}_{scan_cnt_test:03d}.nii.gz')

            test_images.append(str(image_file_nnunet))
            test_labels.append(str(label_file_nnunet))

            # copy the files to new structure
            shutil.copyfile(image_file, image_file_nnunet)

            conversion_dict[str(os.path.abspath(image_file))] = image_file_nnunet

    #Display of number of training and number of testing images
    print("Number of images for training: " + str(scan_cnt_train))
    print("Number of images for testing: " + str(scan_cnt_test))

    #----------------- CREATION OF THE DICTIONNARY-----------------------------------
    # create dataset_description.json
    json_object = json.dumps(conversion_dict, indent=4)
    # write to dataset description
    conversion_dict_name = f"conversion_dict.json"
    with open(os.path.join(path_out, conversion_dict_name), "w") as outfile:
        outfile.write(json_object)


    # c.f. dataset json generation. This contains the metadata for the dataset that nnUNet uses during preprocessing and training
    # general info : https://github.com/MIC-DKFZ/nnUNet/blob/master/nnunet/dataset_conversion/utils.py
    # example: https://github.com/MIC-DKFZ/nnUNet/blob/master/nnunet/dataset_conversion/Task055_SegTHOR.py

    json_dict = OrderedDict()
    json_dict['name'] = args.taskname
    json_dict['description'] = args.taskname
    json_dict['tensorImageSize'] = "3D"
    json_dict['reference'] = "TBD"
    json_dict['licence'] = "TBD"
    json_dict['release'] = "0.0"
    
    # Because only using one modality  
    ## was changed from 'modality' to 'channel_names'
    json_dict['channel_names'] = {
            "0": "T1w",
        }
    
    # 0 is always the background. Any class labels should start from 1.
    json_dict['labels'] = {
        "background" : "0",
        "MS Lesion" : "1" ,
    }
   
    json_dict['numTraining'] = scan_cnt_train
    json_dict['numTest'] = scan_cnt_test
    #Newly required field in the json file with v2
    json_dict["file_ending"] = ".nii.gz"

    json_dict['training'] = [{'image': str(train_labels[i]).replace("labelsTr", "imagesTr") , "label": train_labels[i] }
                                 for i in range(len(train_images))]
    # Note: See https://github.com/MIC-DKFZ/nnUNet/issues/407 for how this should be described

    #Removed because useless in this case
    json_dict['test'] = [str(test_labels[i]).replace("labelsTs", "imagesTs") for i in range(len(test_images))]

    # create dataset_description.json
    json_object = json.dumps(json_dict, indent=4)
    # write to dataset description
    # nn-unet requires it to be "dataset.json"
    dataset_dict_name = f"dataset.json"
    with open(os.path.join(path_out, dataset_dict_name), "w") as outfile:
        outfile.write(json_object)

    return None


if __name__ == '__main__':
    main()