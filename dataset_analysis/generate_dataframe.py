"""
This file generates the dataframe used for the analysis of the CanProCo dataset.
For each patients in the CanProCo dataset, it gathers information about the patient, the pathology and the lesion(s).


Args:
    --data : path to the CanProCo dataset
    --lesion : path to the lesion segmentation file
    --discs : path to the vertebral levels segmentation file
    --spinal-cord : path to the spinal cord segmentation file
    --timepoint : timepoint of the analysis (M0, M12)
    --exclude-file : a .yml file containing the list of participants to exclude
    --output : path to the output folder

Returns:
    None

Example:
    python generate_dataframe.py --data /path/to/CanProCo --lesion /path/to/lesion/segmentation --discs /path/to/discs/segmentation --spinal-cord /path/to/spinal/cord/segmentation --timepoint M0 --exclude-file /path/to/exclude/file --output /path/to/output/folder

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
import yaml
from image import Image, get_dimension, change_orientation


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
    parser.add_argument('--lesion', '-l', type=str, help='path to the lesion segmentation file')
    parser.add_argument('--discs', type=str, help='path to the vertebral levels segmentation file')
    parser.add_argument('--spinal-cord', type=str, help='path to the spinal cord segmentation file')
    parser.add_argument('--timepoint', '-t', type=str, help='timepoint of the analysis (M0, M12)')
    parser.add_argument('--exclude-file', '-e', type=str, help='a .yml file containing the list of participants to exclude')
    parser.add_argument('--output', '-o', type=str, help='path to the output folder')

    return parser


def get_spinal_cord_info(patient_data, spinal_cord_path, timepoint):
    """
    This functions computes the volume of the spinal cord for each patient.
    The volume is added to the patient_data dictionary.
    
    Input:
        patient_data : dictionary containing information about the patient
        dataset_path : path to the CanProCo dataset
        timepoint : timepoint of the analysis (M0, M12)
    
    Returns:
        patient_data : dictionary containing information about the patient
    """

    #get the participant_id
    participant_id = patient_data["participant_id"]

    #now we find the spinal cord segmentation file
    sc_seg_file = os.path.join(spinal_cord_path, participant_id, "ses-" + timepoint, "anat", f"{participant_id}_ses-" + timepoint+ "_PSIR_seg-manual.nii.gz")
    if not os.path.exists(sc_seg_file):
        # If PSIR doesn't exist, use STIR
        sc_seg_file = os.path.join(spinal_cord_path, participant_id, "ses-" + timepoint, "anat", f"{participant_id}_ses-" + timepoint+ "_STIR_seg-manual.nii.gz")
    
    #now we read the spinal cord segmentation file
    sc_seg = nib.load(sc_seg_file)

    #now we get the total volume of the spinal cord
    sc_seg_data = sc_seg.get_fdata()
    voxel_size = sc_seg.header.get_zooms()
    sc_volume = np.sum(sc_seg_data)*voxel_size[0]*voxel_size[1]*voxel_size[2]

    #we add this information to the patient_data dictionary
    patient_data["sc_volume"] = sc_volume

    return patient_data


def analyse_lesion_per_levels(patient_data, discs_path, timepoint, output_folder):
    """
    This function focuses on lesions per vertebral levels.
    It computes the number, volume and average length of lesions per vertebral levels.

    Input:
        lesion_seg_file : path to the lesion segmentation file
        labelled_sc_seg_file : path to the labelled spinal cord segmentation file
        timepoint : timepoint of the analysis (M0, M12)
        output_folder : path to the output folder where the analysis will be saved

    Returns:
        patient_data : dictionary containing information about the patient
    """
    #get the participant_id
    participant_id = patient_data["participant_id"]

    #the output folder
    output_folder = os.path.join(output_folder, "tmp", participant_id)

    #now we find the labelled lesion segmentation file (it was created at the previous step with sct_analyze lesions) and the disc levels segmentation file
    lesion_file_PSIR = os.path.join(output_folder, f"{participant_id}_ses-" + timepoint+ "_PSIR_lesion-manual_label.nii.gz")
    if os.path.exists(lesion_file_PSIR):
        #if PSIR exists, use PSIR
        lesion_seg_file = lesion_file_PSIR
        disc_seg_file = os.path.join(discs_path, participant_id, "ses-" + timepoint, "anat", f"{participant_id}_ses-" + timepoint+ "_PSIR_labels-disc.nii.gz")
    else:
        # If PSIR doesn't exist, use STIR
        lesion_seg_file = os.path.join(output_folder, f"{participant_id}_ses-" + timepoint+ "_STIR_lesion-manual_label.nii.gz")
        disc_seg_file = os.path.join(discs_path, participant_id, "ses-" + timepoint, "anat", f"{participant_id}_ses-" + timepoint+ "_STIR_labels-disc.nii.gz")
    
    #now we read the lesion segmentation file and the disc segmentation file
    lesion_labelled_seg = Image(lesion_seg_file)
    
    #we check if the lesion segmentation file is empty
    if np.sum(lesion_labelled_seg.data) == 0:
        return patient_data

    disc_seg_file = Image(disc_seg_file)
    #we make sure that they have the same orientation
    disc_seg_file = change_orientation(disc_seg_file, lesion_labelled_seg.orientation)

    #we get the data from the lesion segmentation file and the disc segmentation file
    lesion_labelled_seg_data = lesion_labelled_seg.data
    disc_seg_data = disc_seg_file.data

    #get voxel size
    _,_,_,_,voxel_size_x,voxel_size_y,voxel_size_z,_ = get_dimension(lesion_labelled_seg)
    _,_,_,_,disc_voxel_size_x,disc_voxel_size_y,disc_voxel_size_z,_ = get_dimension(disc_seg_file)

    #the levels are the integer unique values without 0
    levels = np.unique(disc_seg_data)
    levels = levels[levels != 0]

    #we iterate over the levels
    for i in range(len(levels)-1):
        #get upper bound of level
        upper_bound = np.where(disc_seg_data == levels[i])
        upper_bound = int(upper_bound[1])*disc_voxel_size_y
        #get lower bound of level
        lower_bound = np.where(disc_seg_data == levels[i+1])
        lower_bound = int(lower_bound[1])*disc_voxel_size_y

        #initialize the lesion_info_per_levels dictionary
        nb_lesions_in_level = 0
        total_lesion_volume_in_level = 0
        lesion_length_in_level = 0

        #now we look for the lesions in this level
        ## get the lesions
        lesions = np.unique(lesion_labelled_seg_data)
        lesions = lesions[lesions != 0]
        ## iterate over the lesions
        for lesion in lesions:
            lesion_voxel = np.where(lesion_labelled_seg_data == lesion)
            center_of_lesion = np.mean(lesion_voxel, axis=1)
            volume_of_lesion = len(lesion_voxel[0])*voxel_size_x*voxel_size_y*voxel_size_z
            top_of_lesion = np.max(lesion_voxel[1])
            bottom_of_lesion = np.min(lesion_voxel[1])
            length_of_lesion = (top_of_lesion - bottom_of_lesion)*voxel_size_y
            #if lesion center is between upper and lower bound, then it is in the level
            if center_of_lesion[1]*voxel_size_y <= upper_bound and center_of_lesion[1]*voxel_size_y >= lower_bound:
                nb_lesions_in_level += 1
                total_lesion_volume_in_level += volume_of_lesion
                lesion_length_in_level += length_of_lesion

        #compute the average vertical length and volume
        avg_lesion_length_in_level = 0
        avg_lesion_volume_in_level = 0
        if nb_lesions_in_level != 0:
            avg_lesion_length_in_level = lesion_length_in_level/nb_lesions_in_level
            avg_lesion_volume_in_level = total_lesion_volume_in_level/nb_lesions_in_level

        #we add this information to the patient_data dictionary
        patient_data[f"nb_lesions_between_{int(levels[i]):02d}_and_{int(levels[i+1]):02d}"] = nb_lesions_in_level
        patient_data[f"total_lesion_volume_between_{int(levels[i]):02d}_and_{int(levels[i+1]):02d}"] = total_lesion_volume_in_level
        patient_data[f"avg_lesion_volume_between_{int(levels[i]):02d}_and_{int(levels[i+1]):02d}"] = avg_lesion_volume_in_level
        patient_data[f"avg_lesion_length_between_{int(levels[i]):02d}_and_{int(levels[i+1]):02d}"] = avg_lesion_length_in_level

        #we add the lesions which are above the 1st level
        if int(levels[i]) == 1:
            #we only have a lower bound
            lower_bound = np.where(disc_seg_data == levels[i])
            lower_bound = int(lower_bound[1])*disc_voxel_size_y

            #initialize the lesion_info_per_levels dictionary
            nb_lesions_above_1 = 0
            total_lesion_volume_above_1 = 0
            lesion_length_above_1 = 0

            #now we look for the lesions in this level
            ## get the lesions
            lesions = np.unique(lesion_labelled_seg_data)
            lesions = lesions[lesions != 0]
            ## iterate over the lesions
            for lesion in lesions:
                lesion_voxel = np.where(lesion_labelled_seg_data == lesion)
                center_of_lesion = np.mean(lesion_voxel, axis=1)
                volume_of_lesion = len(lesion_voxel[0])*voxel_size_x*voxel_size_y*voxel_size_z
                top_of_lesion = np.max(lesion_voxel[1])
                bottom_of_lesion = np.min(lesion_voxel[1])
                length_of_lesion = (top_of_lesion - bottom_of_lesion)*voxel_size_y
                #if lesion center is between upper and lower bound, then it is in the level
                if center_of_lesion[1]*voxel_size_y >= lower_bound:
                    nb_lesions_above_1 += 1
                    total_lesion_volume_above_1 += volume_of_lesion
                    lesion_length_above_1 += length_of_lesion

            #compute the average vertical length and volume
            avg_lesion_length_above_1 = 0
            avg_lesion_volume_above_1 = 0
            if nb_lesions_above_1 != 0:
                avg_lesion_length_above_1 = lesion_length_above_1/nb_lesions_above_1
                avg_lesion_volume_above_1 = total_lesion_volume_above_1/nb_lesions_above_1

            # add info to dictionary
            patient_data[f"nb_lesions_above_01"] = nb_lesions_above_1
            patient_data[f"total_lesion_volume_above_01"] = total_lesion_volume_above_1
            patient_data[f"avg_lesion_volume_above_01"] = avg_lesion_volume_above_1
            patient_data[f"avg_lesion_length_above_01"] = avg_lesion_length_above_1

    # we now add the lesions which are below the last level (we consider that their between the last level and last level + 1)
    i=len(levels)-1
    level = levels[i]
    # we only have an upper bound
    upper_bound = np.where(disc_seg_data == level)
    upper_bound = int(upper_bound[1])*disc_voxel_size_y

    #initialize the lesion_info_per_levels dictionary
    nb_lesions_after_last = 0
    total_lesion_volume_after_last = 0
    lesion_length_after_last = 0

    #now we look for the lesions in this level
    ## get the lesions
    lesions = np.unique(lesion_labelled_seg_data)
    lesions = lesions[lesions != 0]
    ## iterate over the lesions
    for lesion in lesions:
        lesion_voxel = np.where(lesion_labelled_seg_data == lesion)
        center_of_lesion = np.mean(lesion_voxel, axis=1)
        volume_of_lesion = len(lesion_voxel[0])*voxel_size_x*voxel_size_y*voxel_size_z
        top_of_lesion = np.max(lesion_voxel[1])
        bottom_of_lesion = np.min(lesion_voxel[1])
        length_of_lesion = (top_of_lesion - bottom_of_lesion)*voxel_size_y
        #if lesion center is between upper and lower bound, then it is in the level
        if center_of_lesion[1]*voxel_size_y <= upper_bound:
            nb_lesions_after_last += 1
            total_lesion_volume_after_last += volume_of_lesion
            lesion_length_after_last += length_of_lesion

    #compute the average vertical length and volume
    avg_lesion_length_after_last = 0
    avg_lesion_volume_after_last = 0
    if nb_lesions_after_last != 0:
        avg_lesion_length_after_last = lesion_length_after_last/nb_lesions_after_last
        avg_lesion_volume_after_last = total_lesion_volume_after_last/nb_lesions_after_last

    # add info to dictionary
    patient_data[f"nb_lesions_between_{int(levels[i]):02d}_and_{int(levels[i])+1:02d}"] = nb_lesions_after_last
    patient_data[f"total_lesion_between_{int(levels[i]):02d}_and_{int(levels[i])+1:02d}"] = total_lesion_volume_after_last
    patient_data[f"avg_lesion_between_{int(levels[i]):02d}_and_{int(levels[i])+1:02d}"] = avg_lesion_volume_after_last
    patient_data[f"avg_lesion_between_{int(levels[i]):02d}_and_{int(levels[i])+1:02d}"] = avg_lesion_length_after_last

    return patient_data


def analyze_patient_lesion(patient_data, lesion_path, timepoint, output_folder):
    """
    This function analyzes the lesions of a patient using the sct_analyse_lesions.
    It gathers information on the number of lesions, their volumes, their lengths and their equivalent diameters.

    Input:
        patient_data : dictionary containing information about the patient
        dataset_path : path to the CanProCo dataset
        timepoint : timepoint of the analysis (M0, M12)
        output_folder : path to the output folder where the analysis will be saved
    
    Returns:
        patient_data : dictionary containing information about the patient
    """
    #get the participant_id
    participant_id = patient_data["participant_id"]

    #now we find the lesion segmentation file and the spinal cord segmentation file
    lesion_file_PISR = os.path.join(lesion_path, participant_id, "ses-" + timepoint, "anat", f"{participant_id}_ses-" + timepoint+ "_PSIR_lesion-manual.nii.gz")
    if os.path.exists(lesion_file_PISR):
        lesion_seg_file = lesion_file_PISR
        sc_seg_file = os.path.join(lesion_path, participant_id, "ses-" + timepoint, "anat", f"{participant_id}_ses-" + timepoint + "_PSIR_seg-manual.nii.gz")
    else:
        # If PSIR doesn't exist, use STIR
        lesion_seg_file = os.path.join(lesion_path, participant_id, "ses-" + timepoint, "anat", f"{participant_id}_ses-" + timepoint + "_STIR_lesion-manual.nii.gz")
        sc_seg_file = os.path.join(lesion_path, participant_id, "ses-" + timepoint, "anat", f"{participant_id}_ses-" + timepoint + "_STIR_seg-manual.nii.gz")

    #let's use sct to analyze the lesion
    output_folder = os.path.join(output_folder, "tmp", participant_id)
    os.system(f'sct_analyze_lesion -m {lesion_seg_file} -s {sc_seg_file} -ofolder {output_folder} -v 0')

    #now we read the output file
    output_file = str(lesion_seg_file).split('/')[-1].replace('.nii.gz', '_analysis.xls')
    output_file = os.path.join(output_folder, output_file)
    data = pd.read_excel(output_file)
    
    #we check if the data is empty
    if data.empty:
        #if it is empty we return the patient_data dictionary
        return patient_data

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


def analyze_patient_tsv(participant_id, participants_tsv, timepoint):
    """
    This function gathers information from the participants.tsv file.
    It gathers information on each participant, their pathology (and the material used for image acquisition).

    Input:
        participant_id : id of the participant
        participants_tsv : pandas dataframe containing information about the participants
        timepoint : timepoint of the analysis (M0, M12)
    
    Returns:
        patient_data : dictionary containing information about the patient
    """
    patient_data = {}
    patient_data["participant_id"] = participant_id
    #site
    patient_data["site"] = participant_id.split('-')[1][:3]
    #sex
    patient_data["sex"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["sex"].values[0]
    #age at M0
    patient_data["age"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["age_" + timepoint].values[0]
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

    # Get the parsed arguments
    data_path = args.data
    lesion_path = args.lesion
    discs_path = args.discs
    spinal_cord_path = args.spinal_cord
    timepoint = args.timepoint # supposed to be written like M0 or M12

    #exclude some participants (list is like this [sub-mon118, sub-mon006])
    exclude_list = []
    with open(args.exclude_file, 'r') as file:
        exclude_list = yaml.load(file, Loader=yaml.FullLoader)
    #only keep the participant_id
    exclude_list = [participant.split('_')[0] for participant in exclude_list]

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

        # We skip the participant if it is in the exclude list
        if participant in exclude_list:
            print(f"Skipping {participant} because it is in the exclude list")
            continue

        #analyze the patient tsv file
        patient_data = analyze_patient_tsv(participant, participants_tsv, timepoint)

        #analyze the patient lesion
        patient_data = analyze_patient_lesion(patient_data, lesion_path, timepoint, output_folder)

        #analyze the patient lesion distributions per levels
        patient_data = analyse_lesion_per_levels(patient_data, discs_path, timepoint, output_folder)

        #analyze the patient spinal cord volume
        patient_data = get_spinal_cord_info(patient_data, spinal_cord_path, timepoint)

        #add the patient to the dataset
        dataframe = pd.concat([dataframe, pd.DataFrame([patient_data])], ignore_index=True)

    #save the dataset in the output folder
    dataframe.to_csv(os.path.join(output_folder, 'dataframe.csv'), index=False)
    print("Dataframe saved in " + os.path.join(output_folder, 'dataframe.csv'))
    return None


if __name__ == '__main__':
    main()