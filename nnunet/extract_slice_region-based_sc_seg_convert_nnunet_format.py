"""
This script is used to extract slices from the 3D images and labels and save them as 2D images and labels.
For each of these slices, the script segments the spinal cord and saves the images and multi-labels in the nnU-Net format.
The images are stored in the nnU-Net format, i.e. in the folders imagesTr and imagesTs.

Input:
    -path-data: path to the data folder
    -label_folder: name of the folder containing the labels
    -path-out: path to the output folder
    -taskname: name of the task
    -tasknumber: number of the task
    -contrasts: contrasts used for the images
    -qc-folder: path to the quality control folder

Output:
    None

To do:
    -

Example:
    python nnunet/extract_slice_region-based_sc_seg_convert_nnunet_format.py --path-data /path/to/data_extracted --path-out /path/to/nnUNet_raw --taskname TASK-NAME --tasknumber DATASET-ID --label-folder labels --contrasts PSIR --qc_folder /path/to/qc/folder  

Pierre-Louis Benveniste
"""

import argparse
import pathlib
from pathlib import Path
import json
import os
import nibabel as nib
import shutil
from collections import OrderedDict
import numpy as np


def get_parser():
    """
    This function parses the command line arguments and returns an argparse object.

    Input:
        None

    Returns:
        parser : argparse object
    """
    parser = argparse.ArgumentParser(description='Extract labeled slices from 3D images and labels.')
    parser.add_argument('--path-data', required=True,
                        help='Path to BIDS structured dataset. Accepts both cross-sectional and longitudinal datasets')
    parser.add_argument('--label-folder', help='Path to the label folders in derivatives', default='labels', type=str)
    parser.add_argument('--path-out', help='Path to output directory.', required=True)
    parser.add_argument('--taskname', default='MSSpineLesion', type=str,
                        help='Specify the task name - usually the anatomy to be segmented, e.g. Hippocampus',)
    parser.add_argument('--tasknumber', default=501,type=int, 
                        help='Specify the task number, has to be greater than 500 but less than 999. e.g 502')
    parser.add_argument('--contrasts', help='Contrast used for the images', default='PSIR', type=str)
    parser.add_argument('--qc-folder', help='Path to the quality control folder', required=True)
    return parser


def main():
    """
    This function is the main function of the script. It extracts the annotated slices and converts the data from BIDS to nnU-Net format.

    Input:
        None

    Returns:
        None
    """
    #parameters
    label_threshold = 0.0001

    # parse command line arguments
    parser = get_parser()
    args = parser.parse_args()

    path_in_images = Path(args.path_data)
    label_folder = args.label_folder
    path_in_labels = Path(os.path.join(args.path_data, 'derivatives', label_folder))
    path_out = Path(os.path.join(os.path.abspath(args.path_out), f'Dataset{args.tasknumber}_{args.taskname}'))
    contrasts = args.contrasts.split(',')
    qc_folder = args.qc_folder

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
    

    #Initialise the number of scans in train and in test folder
    scan_cnt_train, scan_cnt_test = 0, 0

    valid_train_imgs = []
    valid_test_imgs = []

    #The image folders
    image_files = []
    for contrast in contrasts:
        image_files += list(path_in_images.rglob(f'*_{contrast}.nii.gz'))
    image_files = sorted(image_files)
    
    for image_file in image_files:

        file_id = str(image_file).rsplit('_',1)[0].split('/')[-1] + "_"
        
        if (file_id in str(labelled_imgs)):
            label_file = [k for k in labelled_imgs if file_id in k][0]

            #open the image and the label
            image = nib.load(image_file)
            image_data = image.get_fdata()
            label = nib.load(label_file)
            label_data = label.get_fdata()
            nb_slices = image_data.shape[2]

            if (np.sum(label_data)!=0 ):

                #we segment the spinal cord
                sc_seg_file_path = os.path.join(qc_folder,f'{args.taskname}_{scan_cnt_train:03d}_sc_seg.nii.gz')
                os.system(f'sct_propseg -i {image_file} -o {sc_seg_file_path} -c t1 -v 0')
                #generate the quality control
                os.system(f'sct_qc -i {image_file} -s {sc_seg_file_path} -d {sc_seg_file_path} -p sct_deepseg_lesion -plane sagittal -qc {qc_folder}')
                
                for slice in range(nb_slices):
                    #we check if the label slice is empty or not
                    label_slice = np.asarray(label.dataobj)[:,:,slice]
                    if np.sum(label_slice)!=0 :
                        #we add the slice to the train folder
                        scan_cnt_train+= 1
                        print("Number of images for training: ", scan_cnt_train)

                        image_file_nnunet = os.path.join(path_out_imagesTr,f'{args.taskname}_{scan_cnt_train:03d}_0000.nii.gz')
                        label_file_nnunet = os.path.join(path_out_labelsTr,f'{args.taskname}_{scan_cnt_train:03d}.nii.gz')
                
                        train_images.append(str(image_file_nnunet))
                        train_labels.append(str(label_file_nnunet))

                        #create the image_slice_file and the label_slice_file
                        image_slice = np.asarray(image.dataobj)[:,:,slice]
                        image_slice = np.expand_dims(image_slice, axis=2)
                        image_slice_file = nib.Nifti1Image(image_slice, affine=image.affine)
                        nib.save(image_slice_file, str(image_file_nnunet))
                        label_slice = np.asarray(label.dataobj)[:,:,slice]
                        label_slice = np.expand_dims(label_slice, axis=2)
                        label_slice = np.where(label_slice > label_threshold, 1, 0)

                        #we create the multi-label
                        sc_seg = nib.load(sc_seg_file_path)
                        sc_seg = np.asarray(sc_seg.dataobj)[:,:,slice]
                        sc_seg = np.where(sc_seg > label_threshold, 1, 0)
                        
                        multi_label_slice = np.zeros(image_slice.shape)
                        multi_label_slice[sc_seg==1] = 1
                        multi_label_slice[label_slice==1] = 2

                        multi_label_slice_file = nib.Nifti1Image(multi_label_slice, affine=label.affine)
                        nib.save(multi_label_slice_file, str(label_file_nnunet))

                        conversion_dict[str(os.path.abspath(image_file))] = image_file_nnunet
                        conversion_dict[str(os.path.abspath(label_file))] = label_file_nnunet
        else:

            #open the image and the label
            image = nib.load(image_file)
            image_data = image.get_fdata()
            nb_slices = image_data.shape[2]

            for slice in range(nb_slices):
                #we add the slice to the test folder
                scan_cnt_test += 1
                # create the new convention names
                image_file_nnunet = os.path.join(path_out_imagesTs,f'{args.taskname}_{scan_cnt_test:03d}_0000.nii.gz')
                label_file_nnunet = os.path.join(path_out_labelsTs,f'{args.taskname}_{scan_cnt_test:03d}.nii.gz')

                test_images.append(str(image_file_nnunet))
                test_labels.append(str(label_file_nnunet))

                #create the image_slice_file
                image_slice = np.asarray(image.dataobj)[:,:,slice]
                image_slice = np.expand_dims(image_slice, axis=2)
                image_slice_file = nib.Nifti1Image(image_slice, affine=image.affine)
                nib.save(image_slice_file, str(image_file_nnunet))

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
    json_dict['tensorImageSize'] = "2D"
    json_dict['reference'] = "TBD"
    json_dict['licence'] = "TBD"
    json_dict['release'] = "0.0"
    
    # Because only using one modality  
    ## was changed from 'modality' to 'channel_names'
    json_dict['channel_names'] = {
            "0": "PSIR",
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