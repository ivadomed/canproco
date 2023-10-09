"""
Convert data from BIDS to nnU-Net format
This python script converts data from the BIDS format to the nnU-Net format in order to be able to perform pre-processing, training and inference.
This script prepares the data for the training of a region-based nnunet model predicting MS lesions and spinal cord segmentation.
This script was adapted to build a test set with a given ratio of each site.

Example of run:
    $ python convert_bids_to_nnunet.py --path-data /path/to/data_extracted --path-out /path/to/nnUNet_raw --taskname TASK-NAME --tasknumber DATASET-ID  --contrasts PSIR,STIR --test-ratio 0.2 --include-unlabelled

Arguments:
    --path-data : Path to BIDS structured dataset. Accepts both cross-sectional and longitudinal datasets
    --path-out : Path to output directory.
    --taskname: Specify the task name - usually the anatomy to be segmented, e.g. Hippocampus
    --tasknumber : Specify the task number, has to be greater than 100 but less than 999
    --contrasts : Contrasts used for the images (default='PSIR') (separated by a comma)
    --test-ratio : Ratio of the data to be used for testing (default=0.2)
    --include-unlabelled : Include unlabelled images in the test set (default=False)


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
    parser.add_argument('--contrasts', help='Contrast used for the images', default='PSIR', type=str)
    parser.add_argument('--test-ratio', help='Ratio of the data to be used for testing', default=0.2, type=float)
    parser.add_argument('--include-unlabelled', help='Include unlabelled images in the test set', default=False, type=bool)
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
    path_in_labels = Path(os.path.join(args.path_data, 'derivatives', 'labels'))
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
    
    # We first extract all the label files' names for every wanted contrasask
    label_files = []
    for contrast in contrasts:
        label_files += list(path_in_labels.rglob(f'*_{contrast}_lesion-manual.nii.gz'))
    label_files = sorted(label_files)   
    labelled_imgs += [str(k) for k in label_files]

    ##Number of labelled images per sites
    nb_mon = len([k for k in labelled_imgs if 'sub-mon' in str(k)])
    nb_edm = len([k for k in labelled_imgs if 'sub-edm' in str(k)])
    nb_cal = len([k for k in labelled_imgs if 'sub-cal' in str(k)])
    nb_van = len([k for k in labelled_imgs if 'sub-van' in str(k)])
    nb_tor = len([k for k in labelled_imgs if 'sub-tor' in str(k)])

    #--------------- DISPACTH OF LABELLED IMAGES AND UNLABELED IMAGES ------------------- 
    
    #Initialise the number of scans in train and in test folder
    scan_cnt_train, scan_cnt_test = 0, 0
    #initialise the number of images seen in each site for training
    cnt_mon, cnt_edm, cnt_cal, cnt_van, cnt_tor = 0, 0, 0, 0, 0

    valid_train_imgs = []
    valid_test_imgs = []

    #The image folders
    image_files = []
    for contrast in contrasts:
        image_files += list(path_in_images.rglob(f'*_{contrast}.nii.gz'))
    image_files = sorted(image_files)
    
    for image_file in image_files:

        #Identify common data path
        common = os.path.relpath(os.path.dirname(image_file), args.path_data)

        file_id = str(image_file).rsplit('_',1)[0].split('/')[-1] + "_"
        
        #we check if the image is labelled and if the count is not above the limit set by the ratio
        file_in_mon_and_count = ('sub-mon' in str(file_id)) and (cnt_mon < nb_mon - int(args.test_ratio * nb_mon))
        file_in_edm_and_count = ('sub-edm' in str(file_id)) and (cnt_edm < nb_edm - int(args.test_ratio * nb_edm))
        file_in_cal_and_count = ('sub-cal' in str(file_id)) and (cnt_cal < nb_cal - int(args.test_ratio * nb_cal))
        file_in_van_and_count = ('sub-van' in str(file_id)) and (cnt_van < nb_van - int(args.test_ratio * nb_van))
        file_in_tor_and_count = ('sub-tor' in str(file_id)) and (cnt_tor < nb_tor - int(args.test_ratio * nb_tor))
                                 

        if (file_id in str(labelled_imgs)) and ( file_in_cal_and_count or file_in_edm_and_count or file_in_mon_and_count or file_in_tor_and_count or file_in_van_and_count): 

            lesion_mask_file = [k for k in labelled_imgs if file_id in k][0]

            scan_cnt_train+= 1

            image_file_nnunet = os.path.join(path_out_imagesTr,f'{args.taskname}_{scan_cnt_train:03d}_0000.nii.gz')
            label_file_nnunet = os.path.join(path_out_labelsTr,f'{args.taskname}_{scan_cnt_train:03d}.nii.gz')
            
            train_images.append(str(image_file_nnunet))
            train_labels.append(str(label_file_nnunet))

            #we find the sc_mask
            sc_seg_file = str(lesion_mask_file).replace('lesion-manual', 'sc_seg')
            if not os.path.exists(sc_seg_file):
                #return an error if the sc_seg file does not exist
                raise ValueError(f'The spinal cord segmentation file {sc_seg_file} does not exist. Please run the spinal cord segmentation script first.')

            #we create the multilabelled mask 
            label_threshold = 0.5
            lesion_mask = nib.load(lesion_mask_file)
            lesion_affine = lesion_mask.affine
            lesion_mask = np.asarray(lesion_mask.dataobj)
            lesion_mask = np.where(lesion_mask > label_threshold, 1, 0)

            sc_seg = nib.load(sc_seg_file)
            sc_seg = np.asarray(sc_seg.dataobj)
            sc_seg = np.where(sc_seg > label_threshold, 1, 0)
            
            #we create the multi-label
            multi_label = np.zeros(lesion_mask.shape)
            multi_label[sc_seg==1] = 1
            multi_label[lesion_mask==1] = 2

            multi_label_file = nib.Nifti1Image(multi_label, affine=lesion_affine)
            nib.save(multi_label_file, str(label_file_nnunet))
            
            # copy the files to new structure
            shutil.copyfile(image_file, image_file_nnunet)
           
            conversion_dict[str(os.path.abspath(image_file))] = image_file_nnunet
            conversion_dict[str(os.path.abspath(lesion_mask_file))] = label_file_nnunet

            if file_in_cal_and_count:
                cnt_cal += 1
            elif file_in_edm_and_count:
                cnt_edm += 1
            elif file_in_mon_and_count:
                cnt_mon += 1
            elif file_in_tor_and_count:
                cnt_tor += 1
            elif file_in_van_and_count:
                cnt_van += 1
    
        else:
            #if we only want labelled images in the test set
            if not args.include_unlabelled:

                #if the image is labelled:
                if file_id in str(labelled_imgs):

                    lesion_mask_file = [k for k in labelled_imgs if file_id in k][0]

                    scan_cnt_test+= 1

                    image_file_nnunet = os.path.join(path_out_imagesTs,f'{args.taskname}_{scan_cnt_test:03d}_0000.nii.gz')
                    label_file_nnunet = os.path.join(path_out_labelsTs,f'{args.taskname}_{scan_cnt_test:03d}.nii.gz')
                    
                    test_images.append(str(image_file_nnunet))
                    test_labels.append(str(label_file_nnunet))

                    #we find the sc_mask
                    sc_seg_file = str(lesion_mask_file).replace('lesion-manual', 'sc_seg')
                    if not os.path.exists(sc_seg_file):
                        #return an error if the sc_seg file does not exist
                        raise ValueError(f'The spinal cord segmentation file {sc_seg_file} does not exist. Please run the spinal cord segmentation script first.')

                    #we create the multilabelled mask 
                    label_threshold = 0.5
                    lesion_mask = nib.load(lesion_mask_file)
                    lesion_affine = lesion_mask.affine
                    lesion_mask = np.asarray(lesion_mask.dataobj)
                    lesion_mask = np.where(lesion_mask > label_threshold, 1, 0)

                    sc_seg = nib.load(sc_seg_file)
                    sc_seg = np.asarray(sc_seg.dataobj)
                    sc_seg = np.where(sc_seg > label_threshold, 1, 0)
                    
                    #we create the multi-label
                    multi_label = np.zeros(lesion_mask.shape)
                    multi_label[sc_seg==1] = 1
                    multi_label[lesion_mask==1] = 2

                    multi_label_file = nib.Nifti1Image(multi_label, affine=lesion_affine)
                    nib.save(multi_label_file, str(label_file_nnunet))
                    
                    # copy the files to new structure
                    shutil.copyfile(image_file, image_file_nnunet)
                
                    conversion_dict[str(os.path.abspath(image_file))] = image_file_nnunet
                    conversion_dict[str(os.path.abspath(lesion_mask_file))] = label_file_nnunet
            
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

                #if the image is labelled:
                if file_id in str(labelled_imgs):

                    lesion_mask_file = [k for k in labelled_imgs if file_id in k][0]

                    #we find the sc_mask
                    sc_seg_file = str(lesion_mask_file).replace('lesion-manual', 'sc_seg')
                    if not os.path.exists(sc_seg_file):
                        #return an error if the sc_seg file does not exist
                        raise ValueError(f'The spinal cord segmentation file {sc_seg_file} does not exist. Please run the spinal cord segmentation script first.')

                    #we create the multilabelled mask 
                    label_threshold = 0.5
                    lesion_mask = nib.load(lesion_mask_file)
                    lesion_affine = lesion_mask.affine
                    lesion_mask = np.asarray(lesion_mask.dataobj)
                    lesion_mask = np.where(lesion_mask > label_threshold, 1, 0)

                    sc_seg = nib.load(sc_seg_file)
                    sc_seg = np.asarray(sc_seg.dataobj)
                    sc_seg = np.where(sc_seg > label_threshold, 1, 0)
                    
                    #we create the multi-label
                    multi_label = np.zeros(lesion_mask.shape)
                    multi_label[sc_seg==1] = 1
                    multi_label[lesion_mask==1] = 2

                    multi_label_file = nib.Nifti1Image(multi_label, affine=lesion_affine)
                    nib.save(multi_label_file, str(label_file_nnunet))
                    
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
            "0": "PSIR,STIR",
        }
    
     # 0 is always the background. Any class labels should start from 1.
    json_dict['labels'] = {
        "background" : 0,
        "Spinal cord" : [1, 2] ,
        "Lesion" : [2]
    }
   
    json_dict['regions_class_order'] = [1,2]

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