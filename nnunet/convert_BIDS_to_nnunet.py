"""
Convert data from BIDS to nnU-Net format
This python script converts data from the BIDS format to the nnU-Net format in order to be able to perform pre-processing, training and inference.
This script prepares the data for the training of a region-based nnunet model predicting MS lesions and spinal cord segmentation.
This script was adapted to build a test set with a given ratio of each site.
It can also build the test set with only unlabelled images.

Example of run:
    $ python convert_bids_to_nnunet.py --path-data /path/to/data_extracted --path-out /path/to/nnUNet_raw --taskname TASK-NAME --tasknumber DATASET-ID  --contrasts PSIR,STIR --test-ratio 0.2 --time-point ses-M0 --type training --exclude-file /path/to/exclude_file.yml

Arguments:
    --path-data : Path to BIDS structured dataset. Accepts both cross-sectional and longitudinal datasets
    --path-out : Path to output directory.
    --taskname: Specify the task name - usually the anatomy to be segmented, e.g. Hippocampus
    --tasknumber : Specify the task number, has to be greater than 100 but less than 999
    --contrasts : Contrasts used for the images (default='PSIR') (separated by a comma)
    --test-ratio : Ratio of the data to be used for testing (default=0.2)
    --time-point : Time point of the data to be used (default=ses-M0)
    --type : Type of use of the data (default=training)
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
from tqdm import tqdm

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
    parser.add_argument('--type', help='Type of use of the data (possible answers: "training" or "inference")', default='training', type=str)
    parser.add_argument('--exclude-file', help='Path to the file containing the list of subjects to exclude from the dataset', default=None, type=str)    

    return parser


def  create_multi_label_mask(lesion_mask_file, sc_seg_file, disc_level_file, label_file_nnunet):
    """
    This function creates a multi-label mask for a given lesion mask and spinal cord segmentation mask.
    It also removes the lesions and the spinal cord segmentation which are above the first vertebral level.
    It saves the multi-label mask in the destination folder.

    Input:
        lesion_mask_file : Path to the lesion mask file
        sc_seg_file : Path to the spinal cord segmentation mask file
        label_file_nnunet : Path to the multi-label mask file

    Returns:
        None
    """
    #we create the multilabelled mask 
    label_threshold = 0.5

    #we load the disc levels mask
    disc_level = nib.load(disc_level_file)
    # we get the coordinate of the first disc (label=1)
    first_disc_level = np.where(np.asarray(disc_level.dataobj) == 1)[1]

    #we load the lesion mask
    lesion_mask = nib.load(lesion_mask_file)
    lesion_affine = lesion_mask.affine
    lesion_mask = np.asarray(lesion_mask.dataobj)
    lesion_mask = np.where(lesion_mask > label_threshold, 1, 0)

    #we load the disc levels mask
    sc_seg = nib.load(sc_seg_file)
    sc_seg = np.asarray(sc_seg.dataobj)
    sc_seg = np.where(sc_seg > label_threshold, 1, 0)
    
    #we create the multi-label
    multi_label = np.zeros(lesion_mask.shape)
    multi_label[sc_seg==1] = 1
    multi_label[lesion_mask==1] = 2
    #remove annotations above the first disc level
    multi_label[:,int(first_disc_level) -1:,:] = 0

    #we save it in the destination folder
    multi_label_file = nib.Nifti1Image(multi_label, affine=lesion_affine)
    nib.save(multi_label_file, str(label_file_nnunet))


def build_dataset_for_training(args):
    """
    This functions builds a dataset for training in the nnunet format.
    It builds a dataset with a given ratio of each site, for testing and for training.
    It uses only the contrast and the time point given in the arguments.
    Because we only focus on the spinal cord lesions, we remove the lesion and the sc seg which are above the first vertebral level. 

    Input:
        path_in_images : Path to the images
        path_in_labels : Path to the labels
        path_out : Path to the output directory
        contrasts : List of contrasts
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
    
    cnt_sites = [0 for k in sites]

    valid_train_imgs = []
    valid_test_imgs = []

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

        #we check if the file is labelled 
        if file_id not in str(lesion_masks):
            #then we skip this image
            print("skipping because unlabelled:, ", file_id)
            continue
        
        #we check if the file is in the exclude list
        if file_id.rsplit("_", 1)[0] in exclude_list:
            #then we skip this image
            print("skipping because in exclude list:, ", file_id)
            continue
        
        # find the corresponding label file
        lesion_mask_file = [k for k in lesion_masks if file_id in k][0]

        #we check if the image is not in a site where the number of images is already above the limit ratio (ratio*nb_subj_per_site)
        count_below_limit = count_subj_per_site_dict[file_site] < (1-test_ratio) * nb_subj_per_site_dict[file_site]

        if count_below_limit :
            
            #we update the count of images of the training set
            scan_cnt_train+= 1
            #we update the count of images of the site
            count_subj_per_site_dict[file_site] += 1

            image_file_nnunet = os.path.join(path_out_imagesTr,f'{args.taskname}_{scan_cnt_train:03d}_0000.nii.gz')
            label_file_nnunet = os.path.join(path_out_labelsTr,f'{args.taskname}_{scan_cnt_train:03d}.nii.gz')
            
            train_images.append(str(image_file_nnunet))
            train_labels.append(str(label_file_nnunet))

            #we find the sc mask
            sc_seg_file = str(lesion_mask_file).replace('lesion-manual', 'sc_seg')
            if not os.path.exists(sc_seg_file):
                #return an error if the sc_seg file does not exist
                raise ValueError(f'The spinal cord segmentation file {sc_seg_file} does not exist. Please run the spinal cord segmentation script first.')
            
            #we find the disc levels mask 
            disc_level_file = str(lesion_mask_file).replace('lesion-manual', 'labels-disc')
            if not os.path.exists(disc_level_file):
                #return an error if the disc_level file does not exist
                raise ValueError(f'The disc level file {disc_level_file} does not exist. Please run the disc level script first.')

            #we create the multi-label mask
            create_multi_label_mask(lesion_mask_file, sc_seg_file, disc_level_file, label_file_nnunet)
            
            # copy the image to new structure
            shutil.copyfile(image_file, image_file_nnunet)
           
            #we update the conversion dict (for label we only point to the lesion mask)
            conversion_dict[str(os.path.abspath(image_file))] = image_file_nnunet
            conversion_dict[str(os.path.abspath(lesion_mask_file))] = label_file_nnunet

        else:
            #then we add the image to the test set
            scan_cnt_test += 1
            #we update the count of images of the site
            count_subj_per_site_dict[file_site] += 1

            # create the new convention names
            image_file_nnunet = os.path.join(path_out_imagesTs,f'{args.taskname}_{scan_cnt_test:03d}_0000.nii.gz')
            label_file_nnunet = os.path.join(path_out_labelsTs,f'{args.taskname}_{scan_cnt_test:03d}.nii.gz')
            
            test_images.append(str(image_file_nnunet))
            test_labels.append(str(label_file_nnunet))

            #we find the sc_mask
            sc_seg_file = str(lesion_mask_file).replace('lesion-manual', 'sc_seg')
            if not os.path.exists(sc_seg_file):
                #return an error if the sc_seg file does not exist
                raise ValueError(f'The spinal cord segmentation file {sc_seg_file} does not exist. Please run the spinal cord segmentation script first.')

            #we create the multilabelled masK
            create_multi_label_mask(lesion_mask_file, sc_seg_file, label_file_nnunet)
            
            # copy the files to new structure
            shutil.copyfile(image_file, image_file_nnunet)
        
            conversion_dict[str(os.path.abspath(image_file))] = image_file_nnunet
            conversion_dict[str(os.path.abspath(lesion_mask_file))] = label_file_nnunet
        
        break
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


def build_dataset_for_inference(args):
    """
    This script builds a dataset for inference in the nnunet format.
    It uses only the contrast and the time point given in the arguments.
    It only builds one folder with all the images in the correct naming convention.
    """

    #------------- DEFINITION OF THE PATHS --------------------------
    path_in_images = Path(args.path_data)
    path_out = Path(os.path.join(os.path.abspath(args.path_out), f'Dataset{args.tasknumber}_{args.taskname}'))
    contrasts = args.contrasts.split(',')
    time_point = args.time_point

    #get the list of subjects to exclude
    exclude_list = []
    if args.exclude_file is not None:
       with open(args.exclude_file, 'r') as file:
            exclude_list = yaml.load(file, Loader=yaml.FullLoader)

    #create the output folder
    pathlib.Path(path_out).mkdir(parents=True, exist_ok=True)

    #initialise the conversion dict
    conversion_dict = {}

    #We get the list of all images with given time point and contrast
    image_files = []
    for contrast in contrasts:
        image_files += list(path_in_images.rglob(f'*{time_point}_{contrast}.nii.gz'))
    image_files = sorted(image_files)

    #initialise the count of images
    scan_cnt = 0

    #we iterate over all images
    for image_file in image_files:
            
            file_id = str(image_file).split('/')[-1].split('.')[0]
            
            #we check if the file is in the exclude list
            if file_id.rsplit("_", 1)[0] in exclude_list:
                #then we skip this image
                print("skipping because in exclude list:, ", file_id)
                continue
            
            #we update the count of images of the training set
            scan_cnt+= 1
            
            # create the new convention names
            image_file_nnunet = os.path.join(path_out,f'{args.taskname}_{scan_cnt:03d}.nii.gz')
    
            # copy the files to new structure
            shutil.copyfile(image_file, image_file_nnunet)
            
            conversion_dict[str(os.path.abspath(image_file))] = image_file_nnunet

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

    #we check if we want to build a dataset for nnunet training or for inference
    if args.type == 'training':
        #we build a dataset for training in the nnunet format
        build_dataset_for_training(args)
    elif args.type == 'inference':
        #we build a dataset for inference in the nnunet format
        build_dataset_for_inference(args)
    
    return None


if __name__ == '__main__':
    main()