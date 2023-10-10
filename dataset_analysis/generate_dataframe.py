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


def analyze_patient(participant_id, data_path, participants_tsv, output_folder):
    """
    This function analyzes a patient in the CanProCo dataset.
    It gathers information on each participant, their pathology and their lesion(s).

    Input:
        participant_id : id of the participant
        data_path : path to the CanProCo dataset
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
    # #date_of_scan_M0
    # patient_data["date_of_scan"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["date_of_scan_M0"].values[0]
    # #institution_id_M0
    # patient_data["institution_id"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["institution_id_M0"].values[0]
    # #institution_M0
    # patient_data["institution"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["institution_M0"].values[0]
    # #manufacturer_M0
    # patient_data["manufacturer"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["manufacturer_M0"].values[0]
    # #manufacturers_model_name_M0
    # patient_data["manufacturers_model_name"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["manufacturers_model_name_M0"].values[0]
    # #receive_coil_name_M0
    # patient_data["receive_coil_name"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["receive_coil_name_M0"].values[0]
    # #software_versions_M0
    # patient_data["software_versions"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["software_versions_M0"].values[0]
    #edss_M0
    patient_data["edss"] = participants_tsv.loc[participants_tsv["participant_id"] == participant_id]["edss_M0"].values[0]

    #now we find the lesion segmentation file
    lesion_file_PISR = os.path.join(data_path,"derivatives", "labels", participant_id, "ses-M0", "anat", f"{participant_id}_ses-M0_PSIR_lesion-manual.nii.gz")
    if os.path.exists(lesion_file_PISR):
        lesion_seg_file = lesion_file_PISR
    else:
        # If PSIR doesn't exist, use STIR
        lesion_file_STIR = os.path.join(data_path, "derivatives", "labels", participant_id, "ses-M0", "anat", f"{participant_id}_ses-M0_STIR_lesion-manual.nii.gz")
        lesion_seg_file = lesion_file_STIR

    #we find the spinal cord segmentation file
    sc_seg_file_PSIR = os.path.join(data_path, "derivatives", "labels", participant_id, "ses-M0", "anat", f"{participant_id}_ses-M0_PSIR_sc_seg.nii.gz")
    if os.path.exists(sc_seg_file_PSIR):
        sc_seg_file = sc_seg_file_PSIR
    else:
        # If PSIR doesn't exist, use STIR
        sc_seg_file_STIR = os.path.join(data_path, "derivatives", "labels", participant_id, "ses-M0", "anat", f"{participant_id}_ses-M0_STIR_sc_seg.nii.gz")
        sc_seg_file = sc_seg_file_STIR

    #let's use sct to analyze the lesion
    output_folder = os.path.join(output_folder, "tmp")
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


def main():
    """
    This function generates the dataset used for the analysis of the CanProCo dataset.

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
    dataframe = pd.DataFrame(columns=[
        'participant_id', 'sex', 'age', 'pathology', 'phenotype',
        'edss', 'number_of_lesions', 'total_lesion_volume',
        'biggest_lesion_vol', 'biggest_lesion_length', 'biggest_lesion_eq_diam'
    ])

    #iterate over participant_id
    for participant in participants_tsv["participant_id"]:

        #analyze the patient
        patient_data = analyze_patient(participant, data_path, participants_tsv, output_folder)

        #add the patient to the dataset
        dataframe = pd.concat([dataframe, pd.DataFrame([patient_data])], ignore_index=True)

    print(dataframe)
    #save the dataset in the output folder
    dataframe.to_csv(os.path.join(output_folder, 'dataframe.csv'), index=False)



if __name__ == '__main__':
    main()