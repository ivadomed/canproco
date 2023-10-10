"""
This file generates the dataframe used for the analysis of the CanProCo dataset.
For each patients in the CanProCo dataset, it gathers information about the patient, the pathology and the lesion(s).
For now, this only focuses on the M0 timepoint.

Args:
    --data : path to the CanProCo dataset
    --output : path to the output folder

Returns:
    None

Example:
    python generate_dataframe.py --data /path/to/CanProCo --output /path/to/output

Todo:
    *

Pierre-Louis Benveniste
"""

import os
import argparse
from pathlib import Path
import pandas as pd
import json
import shutil
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
    parser = argparse.ArgumentParser(description='Generates the dataset used for the analysis of the CanProCo dataset.')
    parser.add_argument('--data', '-d', type=str, help='path to the CanProCo dataset')
    parser.add_argument('--output', '-o', type=str, help='path to the output folder')

    return parser


def analyse_lesion_per_levels(patient_data, dataset_path, output_folder):
    """
    This function focuses on lesions per spinal cord levels.
    It computes the number, volume and average length of lesions per spinal cord level.

    Input:
        lesion_seg_file : path to the lesion segmentation file
        labelled_sc_seg_file : path to the labelled spinal cord segmentation file

    Returns:
        patient_data : dictionary containing information about the patient
    """
    #get the participant_id
    participant_id = patient_data["participant_id"]

    #the output folder
    output_folder = os.path.join(output_folder, "tmp", participant_id)

    #now we find the labelledlesion segmentation file (it was created at the previous step with sct_analyze lesions) and the vertebral levels segmentation file
    lesion_file_PSIR = os.path.join(output_folder, f"{participant_id}_ses-M0_PSIR_lesion-manual_label.nii.gz")
    if os.path.exists(lesion_file_PSIR):
        #if PSIR exists, use PSIR
        lesion_seg_file = lesion_file_PSIR
        vert_levels_seg_file = os.path.join(dataset_path, "derivatives", "labels", participant_id, "ses-M0", "anat", f"{participant_id}_ses-M0_PSIR_sc_seg_labeled_discs.nii.gz")
    else:
        # If PSIR doesn't exist, use STIR
        lesion_seg_file = os.path.join(output_folder, f"{participant_id}_ses-M0_STIR_lesion-manual_label.nii.gz")
        vert_levels_seg_file = os.path.join(dataset_path, "derivatives", "labels", participant_id, "ses-M0", "anat", f"{participant_id}_ses-M0_STIR_sc_seg_labeled_discs.nii.gz")
    
    #now we read the lesion segmentation file and the labelled spinal cord segmentation file
    lesion_labelled_seg = nib.load(lesion_seg_file)
    vert_levels_seg_file = nib.load(vert_levels_seg_file)

    #we get the data from the lesion segmentation file and the labelled spinal cord segmentation file
    lesion_labelled_seg_data = lesion_labelled_seg.get_fdata()
    vert_levels_seg_data = vert_levels_seg_file.get_fdata()

    #get voxel size
    voxel_size = lesion_labelled_seg.header.get_zooms()

    #the levels are the integer unique values without 0
    levels = np.unique(vert_levels_seg_data)
    levels = levels[levels != 0]

    #we iterate over the levels
    for i in range(len(levels)-1):
        #get upper bound of level
        upper_bound = np.where(vert_levels_seg_data == levels[i])
        upper_bound = int(upper_bound[1])
        #get lower bound of level
        lower_bound = np.where(vert_levels_seg_data == levels[i+1])
        lower_bound = int(lower_bound[1])

        #initialize the lesion_info_per_levels dictionary
        nb_lesions_in_level = 0
        total_lesion_volume_in_level = 0
        avg_lesion_length_in_level = 0

        #now we look for the lesions in this level
        ## get the lesions
        lesions = np.unique(lesion_labelled_seg_data)
        lesions = lesions[lesions != 0]
        ## iterate over the lesions
        for lesion in lesions:
            lesion_voxel = np.where(lesion_labelled_seg_data == lesion)
            center_of_lesion = np.mean(lesion_voxel, axis=1)
            volume_of_lesion = len(lesion_voxel[0])*voxel_size[0]*voxel_size[1]*voxel_size[2]
            top_of_lesion = np.max(lesion_voxel[2])
            bottom_of_lesion = np.min(lesion_voxel[2])
            length_of_lesion = (top_of_lesion - bottom_of_lesion)*voxel_size[2]
            #if lesion center is between upper and lower bound, then it is in the level
            if center_of_lesion[1] <= upper_bound and center_of_lesion[1] >= lower_bound:
                nb_lesions_in_level += 1
                total_lesion_volume_in_level += volume_of_lesion
                avg_lesion_length_in_level += length_of_lesion

        #compute the average vertical length
        if nb_lesions_in_level != 0:
            avg_lesion_length_in_level = avg_lesion_length_in_level/nb_lesions_in_level

        #we add this information to the patient_data dictionary
        patient_data[f"nb_lesions_between_{levels[i]}_and_{levels[i+1]}"] = nb_lesions_in_level
        patient_data[f"total_lesion_volume_between_{levels[i]}_and_{levels[i+1]}"] = total_lesion_volume_in_level
        patient_data[f"avg_lesion_length_between_{levels[i]}_and_{levels[i+1]}"] = avg_lesion_length_in_level

    return patient_data


def analyze_patient_lesion(patient_data, dataset_path, output_folder):
    """
    This function analyzes the lesions of a patient using the sct_analyse_lesions.
    It gathers information on the number of lesions, their volumes, their lengths and their equivalent diameters.

    Input:
        patient_data : dictionary containing information about the patient
        dataset_path : path to the CanProCo dataset
        output_folder : path to the output folder where the analysis will be saved
    
    Returns:
        patient_data : dictionary containing information about the patient
    """
    #get the participant_id
    participant_id = patient_data["participant_id"]

    #now we find the lesion segmentation file and the spinal cord segmentation file
    lesion_file_PISR = os.path.join(dataset_path,"derivatives", "labels", participant_id, "ses-M0", "anat", f"{participant_id}_ses-M0_PSIR_lesion-manual.nii.gz")
    if os.path.exists(lesion_file_PISR):
        lesion_seg_file = lesion_file_PISR
        sc_seg_file = os.path.join(dataset_path, "derivatives", "labels", participant_id, "ses-M0", "anat", f"{participant_id}_ses-M0_PSIR_sc_seg.nii.gz")
    else:
        # If PSIR doesn't exist, use STIR
        lesion_seg_file = os.path.join(dataset_path, "derivatives", "labels", participant_id, "ses-M0", "anat", f"{participant_id}_ses-M0_STIR_lesion-manual.nii.gz")
        sc_seg_file = os.path.join(dataset_path, "derivatives", "labels", participant_id, "ses-M0", "anat", f"{participant_id}_ses-M0_STIR_sc_seg.nii.gz")

    #let's use sct to analyze the lesion
    output_folder = os.path.join(output_folder, "tmp", participant_id)
    os.system(f'sct_analyze_lesion -m {lesion_seg_file} -s {sc_seg_file} -ofolder {output_folder} -v 0')

    #now we read the output file
    output_file = str(lesion_seg_file).split('/')[-1].replace('.nii.gz', '_analysis.xls')
    output_file = os.path.join(output_folder, output_file)
    data = pd.read_excel(output_file)
    number_of_lesions = max(data['label'])
    total_lesion_volume = sum(data['volume [mm3]'])
    biggest_lesion_vol = max(data['volume [mm3]'])
    biggest_lesion_length = max(data['length [mm]'])
    biggest_lesion_eq_diam = max(data['max_equivalent_diameter [mm]'])

    #we add this information to the patient_data dictionary
    patient_data["number_of_lesions"] = number_of_lesions
    patient_data["total_lesion_volume"] = total_lesion_volume
    patient_data["biggest_lesion_vol"] = biggest_lesion_vol
    patient_data["biggest_lesion_length"] = biggest_lesion_length
    patient_data["biggest_lesion_eq_diam"] = biggest_lesion_eq_diam

    return patient_data


def analyze_patient_tsv(participant_id, participants_tsv):
    """
    This function gathers information from the participants.tsv file.
    It gathers information on each participant, their pathology (and the material used for image acquisition).

    Input:
        participant_id : id of the participant
        participants_tsv : pandas dataframe containing information about the participants
    
    Returns:
        patient_data : dictionary containing information about the patient
    """
    patient_data = {}
    patient_data["participant_id"] = participant_id
    #sex
    patient_data["sex"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["sex"].values[0]
    #age at M0
    patient_data["age"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["age_M0"].values[0]
    #pathology
    patient_data["pathology"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["pathology_M0"].values[0]
    #phenotype_M0
    patient_data["phenotype"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["phenotype_M0"].values[0]
    #edss_M0
    patient_data["edss"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["edss_M0"].values[0]

    #### ------------- THE BELOW INFORMATION IS NOT USED FOR NOW -------------- ####
    # # #date_of_scan_M0
    # patient_data["date_of_scan"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["date_of_scan_M0"].values[0]
    # # #institution_id_M0
    # patient_data["institution_id"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["institution_id_M0"].values[0]
    # # #institution_M0
    # patient_data["institution"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["institution_M0"].values[0]
    # # #manufacturer_M0
    # patient_data["manufacturer"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["manufacturer_M0"].values[0]
    # # #manufacturers_model_name_M0
    # patient_data["manufacturers_model_name"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["manufacturers_model_name_M0"].values[0]
    # # #receive_coil_name_M0
    # patient_data["receive_coil_name"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["receive_coil_name_M0"].values[0]
    # # #software_versions_M0
    # patient_data["software_versions"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["software_versions_M0"].values[0]
    #### ------------- ------------------------------------------ -------------- ####
    

    return patient_data


def main():
    """
    This function generates the dataframe used for the analysis of the CanProCo dataset.

    Input:
        None
    
    Returns:
        None
    """

    # Parse the command line arguments
    parser = get_parser()
    args = parser.parse_args()

    # Get the path to the CanProCo dataset
    data_path = args.data

    #build the output folder
    output_folder = args.output
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    #load the participants.tsv file
    participants_file = os.path.join(data_path, 'participants.tsv')
    participants_tsv = pd.read_csv(participants_file, sep='\t')


    # Create an empty DataFrame
    dataframe = pd.DataFrame()

    #iterate over participant_id
    for participant in participants_tsv["participant_id"]:

        #analyze the patient tsv file
        patient_data = analyze_patient_tsv(participant, participants_tsv)

        #analyze the patient lesion
        patient_data = analyze_patient_lesion(patient_data, data_path, output_folder)

        #analyze the patient lesion distributions per levels
        patient_data = analyse_lesion_per_levels(patient_data, data_path, output_folder)

        #add the patient to the dataset
        dataframe = pd.concat([dataframe, pd.DataFrame([patient_data])], ignore_index=True)

    print(dataframe)
    #save the dataset in the output folder
    dataframe.to_csv(os.path.join(output_folder, 'dataframe.csv'), index=False)



if __name__ == '__main__':
    main()