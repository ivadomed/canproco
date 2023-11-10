"""
Convert data from BIDS to nnU-Net format
This python script converts data from the BIDS format to the nnU-Net format in order to be able to perform pre-processing, training and inference.
This script prepares the data for the training of a nnU-Net model taking two channels as input (i) STIR and invertes PSIR, and (ii)  SC mask.
This script was adapted to build a test set with a given ratio of each site.

Example of run:
    python convert_BIDS_to_nnunet_mul_PSIR.py --path-data /path/to/data --path-out /output/path --taskname MSSCLesion --tasknumber 111 --contrasts PSIR,STIR --test-ratio 0.2 --time-point M0 --exclude-file /path/to/exclude_file

Arguments:
    --path-data : Path to BIDS structured dataset. Accepts both cross-sectional and longitudinal datasets
    --path-out : Path to output directory.
    --taskname: Specify the task name - usually the anatomy to be segmented, e.g. Hippocampus
    --tasknumber : Specify the task number, has to be greater than 100 but less than 999
    --contrasts : Contrasts used for the images (default='PSIR') (separated by a comma)
    --test-ratio : Ratio of the data to be used for testing (default=0.2)
    --time-point : Time point of the data to be used (default=ses-M0)
    --exclude-file : Path to the file containing the list of subjects to exclude from the dataset (default=None)


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

import nibabel as nib
import numpy as np
import yaml


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
    parser.add_argument('--contrasts', help='Contrast used for the images', default='PSIR', type=str)
    parser.add_argument('--test-ratio', help='Ratio of the data to be used for testing', default=0.2, type=float)
    parser.add_argument('--time-point', help='Time point of the data to be used', default='M0', type=str)

    return parser


def build_dataset_for_training(args):
    """
    This functions builds a dataset for training in the nnunet format.
    It builds a dataset with a given ratio of each site, for testing and for training.
    It uses only the contrast and the time point given in the arguments.
    Because we only focus on the spinal cord lesions, we remove the lesion and the sc seg which are above the first vertebral level. 

    Input:
        args : Arguments of the script
    
    Returns:
        None
    """
    #------------- DEFINITION OF THE PATHS --------------------------
    path_in_images = Path(args.path_data)
    path_in_labels = Path(os.path.join(args.path_data, 'derivatives', 'labels'))
    path_out = Path(os.path.join(os.path.abspath(args.path_out), f'Dataset{args.tasknumber}_{args.taskname}'))
    contrasts = args.contrasts.split(',')
    time_point = args.time_point
    test_ratio = args.test_ratio

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

    #initialise the conversion dict
    conversion_dict = {}

    #------------- EXTRACTION OF THE LABELED IMAGES NAMES--------------------------
    
    # We first extract all the lesion mask files' names for every wanted contrasask
    lesion_masks = []
    for contrast in contrasts:
        lesion_masks += list(path_in_labels.rglob(f'*{time_point}_{contrast}_lesion-manual.nii.gz'))
    lesion_masks = sorted(lesion_masks)   
    lesion_masks = [str(k) for k in lesion_masks]

    #we then get all the sc seg files' names for every wanted contrast
    sc_masks = []
    for contrast in contrasts:
        sc_masks += list(path_in_labels.rglob(f'*{time_point}_{contrast}_seg-manual.nii.gz'))
    sc_masks = sorted(sc_masks)
    sc_masks = [str(k) for k in sc_masks]

    #get the site names
    subj = [str(k).split('/')[-1].split('_')[0][:-3] for k in lesion_masks]
    sites = np.unique(subj)
    #we count the number of subjects per site
    nb_subj_per_site = [len([k for k in lesion_masks if site in k]) for site in sites]
    #we build a dictionnary with the number of subjects per site
    nb_subj_per_site_dict = dict(zip(sites, nb_subj_per_site))

    #get the list of subjects to exclude
    exclude_list = []
    if args.exclude_file is not None:
       with open(args.exclude_file, 'r') as file:
            exclude_list = yaml.load(file, Loader=yaml.FullLoader)

    #--------------- DISPACTH OF LABELLED IMAGES IN TRAIN OR TEST SET ------------------- 
    
    #Initialise the number of scans in train and in test folder
    scan_cnt_train, scan_cnt_test = 0, 0
    
    #we build a dictionnary with the counts of images per site
    count_subj_per_site_dict = dict(zip(sites, [0 for k in sites]))


    #We get the list of all images with given time point and contrast
    image_files = []
    for contrast in contrasts:
        image_files += list(path_in_images.rglob(f'*{time_point}_{contrast}.nii.gz'))
    image_files = sorted(image_files)
    
    #we iterate over all images
    for image_file in image_files:

        #get file_id (example: sub-mon001_ses-M0_PSIR)
        file_id = str(image_file).split('/')[-1].split('.')[0]
        #we extract the site of the file
        file_site = str(image_file).split('/')[-1].split('_')[0][:-3]

        #we check if the lesion mask exist for this file 
        if file_id not in str(lesion_masks):
            #then we skip this image
            print("skipping because unlabelled:, ", file_id)
            continue

        #we check if the sc mask exist for this file
        if file_id not in str(sc_masks):
            #then we skip this image
            print("skipping because no sc seg:, ", file_id)
            continue
        
        #we check if the file is in the exclude list
        if file_id.rsplit("_", 1)[0] in exclude_list:
            #then we skip this image
            print("skipping because in exclude list:, ", file_id)
            continue
        
        # find the corresponding label file
        lesion_mask_file = [k for k in lesion_masks if file_id in k][0]
        sc_mask_file = [k for k in sc_masks if file_id in k][0]

        #we check if the image is not in a site where the number of images is already above the limit ratio (ratio*nb_subj_per_site)
        count_below_limit = count_subj_per_site_dict[file_site] < (1-test_ratio) * nb_subj_per_site_dict[file_site]

        if count_below_limit :
            
            #we update the count of images of the training set
            scan_cnt_train+= 1
            #we update the count of images of the site
            count_subj_per_site_dict[file_site] += 1

            image_file_nnunet_channel_1 = os.path.join(path_out_imagesTr,f'{args.taskname}_{scan_cnt_train:03d}_0000.nii.gz')
            image_file_nnunet_channel_2 = os.path.join(path_out_imagesTr,f'{args.taskname}_{scan_cnt_train:03d}_0001.nii.gz')
            label_file_nnunet = os.path.join(path_out_labelsTr,f'{args.taskname}_{scan_cnt_train:03d}.nii.gz')
            
            train_images.append(str(image_file_nnunet_channel_1))
            train_images.append(str(image_file_nnunet_channel_2))
            train_labels.append(str(label_file_nnunet))

            # copy the image to new structure
            if 'PSIR' in str(image_file):
                #we use sct_maths to directly saved the image multiplied by -1 in the nnunet dataset
                os.system(f'sct_maths -i {image_file} -mul -1 -o {image_file_nnunet_channel_1}')
                
            else:
                shutil.copyfile(image_file, image_file_nnunet_channel_1)

            #we also copy the sc mask to the nnunet dataset
            shutil.copyfile(sc_mask_file, image_file_nnunet_channel_2)
                
           
            #we update the conversion dict (for label we only point to the lesion mask)
            conversion_dict[str(os.path.abspath(image_file))] = image_file_nnunet_channel_1
            conversion_dict[str(os.path.abspath(sc_mask_file))] = image_file_nnunet_channel_2

            conversion_dict[str(os.path.abspath(lesion_mask_file))] = label_file_nnunet

        else:
            #if not for training, we add the image to the test set
            scan_cnt_test += 1
            #we update the count of images of the site
            count_subj_per_site_dict[file_site] += 1

            # create the new convention names
            image_file_nnunet_channel_1 = os.path.join(path_out_imagesTs,f'{args.taskname}_{scan_cnt_test:03d}_0000.nii.gz')
            image_file_nnunet_channel_2 = os.path.join(path_out_imagesTs,f'{args.taskname}_{scan_cnt_test:03d}_0001.nii.gz')
            label_file_nnunet = os.path.join(path_out_labelsTs,f'{args.taskname}_{scan_cnt_test:03d}.nii.gz')
            
            test_images.append(str(image_file_nnunet_channel_1))
            test_images.append(str(image_file_nnunet_channel_2))
            test_labels.append(str(label_file_nnunet))
            
            # copy the files to new structure
            if 'PSIR' in str(image_file):
                #we use sct_maths to directly saved the image multiplied by -1 in the nnunet dataset
                os.system(f'sct_maths -i {image_file} -mul -1 -o {image_file_nnunet_channel_1}')

            else:
                shutil.copyfile(image_file, image_file_nnunet_channel_1)

            #we also copy the sc mask to the nnunet dataset
            shutil.copyfile(sc_mask_file, image_file_nnunet_channel_2)
        
            conversion_dict[str(os.path.abspath(image_file))] = image_file_nnunet_channel_1
            conversion_dict[str(os.path.abspath(sc_mask_file))] = image_file_nnunet_channel_2
            conversion_dict[str(os.path.abspath(lesion_mask_file))] = label_file_nnunet
        
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
            "0": "PSIR (mul by -1), STIR",
            "1": "SC_seg",
        }
    
     # 0 is always the background. Any class labels should start from 1.
    json_dict['labels'] = {
        "background" : 0,
        "Lesion" : 1,
    }

    json_dict['numTraining'] = scan_cnt_train
    json_dict['numTest'] = scan_cnt_test
    #Newly required field in the json file with v2
    json_dict["file_ending"] = ".nii.gz"

    json_dict['training'] = [{'image': str(train_labels[i]).replace("labelsTr", "imagesTr") , "label": train_labels[i] }
                                 for i in range(len(train_images))]
    # Note: See https://github.com/MIC-DKFZ/nnUNet/issues/407 for how this should be described

    #Removed because useless in this case
    json_dict['test'] = [{'image': str(test_labels[i]).replace("labelsTs", "imagesTs") , "label": test_labels[i] }
                                 for i in range(len(test_images))]

    # create dataset_description.json
    json_object = json.dumps(json_dict, indent=4)
    # write to dataset description
    # nn-unet requires it to be "dataset.json"
    dataset_dict_name = f"dataset.json"
    with open(os.path.join(path_out, dataset_dict_name), "w") as outfile:
        outfile.write(json_object)
    return None


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

    #we build a dataset for training in the nnunet format
    build_dataset_for_training(args)
    
    return None


if __name__ == '__main__':
    main()