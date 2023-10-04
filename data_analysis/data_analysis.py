"""
This python files performs data analysis on the canproco dataset.
It details the number of subjects and the files available for each subject (for different time points and contrasts).
It generates a text file containing the results of the analysis and a json file containing the information for each subject.
It also looks into the participants.tsv file to get the list of healthy controls and see which contrats they have.

Args:
    -d, --dataset-path: path to the dataset
    -o, --output-path: path to the output directory
    -p, --participants-tsv-path: path to the participants.tsv file

Returns:
    - a csv file containing the results of the analysis

Example:
    python data_analysis.py -d /path/to/dataset -o /path/to/output -p /path/to/participants.tsv

To do:
    *

Pierre-Louis Benveniste
"""


import argparse
import os
import json


def get_parser():
    """
    This function parses the arguments given to the script.

    Args:
        None

    Returns:
        parser: parser containing the arguments
    """

    parser = argparse.ArgumentParser(description='Perform data analysis on the canproco dataset')
    parser.add_argument('-d', '--dataset-path', type=str, required=True, help='path to the dataset')
    parser.add_argument('-o', '--output-path', type=str, required=True, help='path to the output directory')
    parser.add_argument('-p', '--participants-tsv-path', type=str, required=True, help='path to the participants.tsv file')
   
    return parser


def get_healthy_control(participants_tsv_path):
    """
    This function gets the list of healthy controls.

    Args:
        participants_tsv_path: path to the participants.tsv file
    
    Returns:
        subjects: list of healthy controls
    """
    #initialize the list of subjects
    subjects = []

    #open the participants.tsv file
    with open(participants_tsv_path, 'r') as f:
        lines = f.readlines()
        #iterate over the lines
        for line in lines[1:]:
            line = line.split('\t')
            #if the subject is a healthy control
            if line[3] == 'HC':
                subjects.append(line[0])
    return subjects


def main():
    """
    This function performs the data analysis.

    Args:
        None
    
    Returns:
        None
    """
    # Get the parser
    parser = get_parser()
    args = parser.parse_args()

    # Get the arguments
    dataset_path = args.dataset_path
    output_path = args.output_path
    participants_tsv_path = args.participants_tsv_path

    # Get the list of subjects
    healthy_controls = get_healthy_control(participants_tsv_path)

    #time points (for now we only work on M0)
    time_points = ['ses-M0', 'ses-M12']

    # Get the list of subjects
    subjects = os.listdir(dataset_path)
    subjects = [subject for subject in subjects if 'sub-' in subject]
    #print("Total number of subjects: {}".format(len(subjects)))

    #initialize lists
    subjects_all_time_points = []
    subjects_no_M0 = []
    subjects_no_M12 = []
    subjects_PSIR = []
    subjects_STIR = []
    subjects_PSIR_STIR = []
    subjects_no_PSIR_no_STIR = []
    subjects_no_PSIR_no_STIR_once = []
    subjects_hc_psir = []
    subjects_hc_stir = []

    subjects_info = {}
    
    #Iterate over the subjects
    for subject in subjects:
        #iterate over the time_points
        #print("Subject: {}".format(subject))
        sub_time_points = []
        for time_point in time_points:
             #if time_point exists for the subject
            if os.path.exists(os.path.join(dataset_path, subject, time_point)):
                sub_time_points.append(time_point)
        #print("Time points available: {}".format(sub_time_points))
        #initialize the contrast_subject dictionary
        contrast_subject = {}
        for time_point in sub_time_points:
            contrast_subject[time_point] = []
        #iterate over the time points
        for time_point in sub_time_points:
            #print("Time point: {}".format(time_point))
            #get the MRI files for the subject
            subject_path = os.path.join(dataset_path, subject, time_point, 'anat')
            subject_files = os.listdir(subject_path)
            subject_files = [file for file in subject_files if '.nii.gz' in file]
            #we get the contrast for each file
            for file in subject_files:
                contrast_subject[time_point].append(file.split('_')[2].split('.')[0])
            #we print the contrasts available for the subject
            #print("Contrasts available: {}".format(sorted(contrast_subject[time_point])))
        #print(contrast_subject)
        #print("-----------------------------------")
        sub_hc = False
        if subject in healthy_controls:
            sub_hc = True
        subject_info = {'subject': subject, 'time_points': sub_time_points, 'contrasts': contrast_subject, 'healthy_control': sub_hc}
        subjects_info[subject] = subject_info

        #we get the list of the subjects with all the time points
        if len(sub_time_points) == len(time_points):
            subjects_all_time_points.append(subject)
        #we get the list of the subjects with no M0
        if 'ses-M0' not in sub_time_points:
            subjects_no_M0.append(subject)
        #we get the list of the subjects with no M12
        if 'ses-M12' not in sub_time_points:
            subjects_no_M12.append(subject)

        #we get the list of the subjects with PSIR at every time point that they have
        psir_present = True
        for time_point in sub_time_points:
            if 'PSIR' not in contrast_subject[time_point]:
                psir_present = False
        if psir_present:
            subjects_PSIR.append(subject)
        #we get the list of the subjects with STIR at every time point that they have
        stir_present = True
        for time_point in sub_time_points:
            if 'STIR' not in contrast_subject[time_point]:
                stir_present = False
        if stir_present:
            subjects_STIR.append(subject)
        #we get the list of the subjects with PSIR and STIR at every time point that they have
        psir_stir_present = True
        for time_point in sub_time_points:
            if 'PSIR' not in contrast_subject[time_point] or 'STIR' not in contrast_subject[time_point]:
                psir_stir_present = False
        if psir_stir_present:
            subjects_PSIR_STIR.append(subject)
        #we get the list of the subjects with no PSIR and no STIR at every time point that they have
        psir_stir_not_present = True
        for time_point in sub_time_points:
            if 'PSIR' in contrast_subject[time_point] or 'STIR' in contrast_subject[time_point]:
                psir_stir_not_present = False
        if psir_stir_not_present:
            subjects_no_PSIR_no_STIR.append(subject)
        #we get the list of the subjects with no PSIR and no STIR at least once
        psir_stir_not_present_once = False
        for time_point in sub_time_points:
            if 'PSIR' not in contrast_subject[time_point] and 'STIR' not in contrast_subject[time_point]:
                psir_stir_not_present_once = True
        if psir_stir_not_present_once:
            subjects_no_PSIR_no_STIR_once.append(subject)
        #we get the number of subjects with PSIR at every time point that they have
        hc_psir_present = True
        for time_point in sub_time_points:
            if 'PSIR' not in contrast_subject[time_point] or subject not in healthy_controls:
                hc_psir_present = False
        if hc_psir_present:
            subjects_hc_psir.append(subject)
        #we get the number of subjects with STIR at every time point that they have
        hc_stir_present = True
        for time_point in sub_time_points:
            if 'STIR' not in contrast_subject[time_point] or subject not in healthy_controls:
                hc_stir_present = False
        if hc_stir_present:
            subjects_hc_psir.append(subject)
        
    #we print the results
    print("Total number of subjects: {}".format(len(subjects)))
    print("Number of healthy controls: {}".format(len(healthy_controls)))
    print("Number of subjects with all time points: {}".format(len(subjects_all_time_points)))
    print("Number of subjects with no M0: {}".format(len(subjects_no_M0)))
    print("Number of subjects with no M12: {}".format(len(subjects_no_M12)))
    print("Number of subjects with PSIR at every time point they have: {}".format(len(subjects_PSIR)))
    print("Number of subjects with STIR at every time point they have: {}".format(len(subjects_STIR)))
    print("Number of subjects with PSIR and STIR at every time point they have: {}".format(len(subjects_PSIR_STIR)))
    print("Number of subjects with no PSIR and no STIR at every time point they have: {}".format(len(subjects_no_PSIR_no_STIR)))
    print("Number of subjects with no PSIR and no STIR at least once: {}".format(len(subjects_no_PSIR_no_STIR_once)))
    print("Number of healthy control with PSIR at every time point they have: {}".format(len(subjects_hc_psir)))
    print("Number of healthy control with STIR at every time point they have: {}".format(len(subjects_hc_stir)))
    print("-----------------------------------")

    #we save the subjects_info dictionary in a json file
    with open(os.path.join(output_path, 'subjects_info.json'), 'w') as fp:
        json.dump(subjects_info, fp, indent=4)
    
    #we write a txt file with the results
    with open(os.path.join(output_path, 'results.txt'), 'w') as f:
        f.write("Total number of subjects: {}\n".format(len(subjects)))
        f.write("Number of healthy controls: {}".format(len(healthy_controls)))
        f.write("Number of subjects with all time points: {}\n".format(len(subjects_all_time_points)))
        f.write("Number of subjects with no M0: {}\n".format(len(subjects_no_M0)))
        f.write("Number of subjects with no M12: {}\n".format(len(subjects_no_M12)))
        f.write("Number of subjects with PSIR at every time point they have: {}\n".format(len(subjects_PSIR)))
        f.write("Number of subjects with STIR at every time point they have: {}\n".format(len(subjects_STIR)))
        f.write("Number of subjects with PSIR and STIR at every time point they have: {}\n".format(len(subjects_PSIR_STIR)))
        f.write("Number of subjects with no PSIR and no STIR at every time point they have: {}\n".format(len(subjects_no_PSIR_no_STIR)))
        f.write("Number of subjects with no PSIR and no STIR at least once: {}\n".format(len(subjects_no_PSIR_no_STIR_once)))
        f.write("Number of healthy control with PSIR at every time point they have: {}\n".format(len(subjects_hc_psir)))
        f.write("Number of healthy control with STIR at every time point they have: {}\n".format(len(subjects_hc_stir)))
        f.write("-----------------------------------\n")
        f.write("Subjects with all time points:\n")
        for subject in subjects_all_time_points:
            f.write("{}\n".format(subject))
        f.write("-----------------------------------\n")
        f.write("Healthy controls:\n")
        for subject in healthy_controls:
            f.write("{}\n".format(subject))
        f.write("-----------------------------------\n")
        f.write("Subjects with no M0:\n")
        for subject in subjects_no_M0:
            f.write("{}\n".format(subject))
        f.write("-----------------------------------\n")
        f.write("Subjects with no M12:\n")
        for subject in subjects_no_M12:
            f.write("{}\n".format(subject))
        f.write("-----------------------------------\n")
        f.write("Subjects with PSIR at every time point they have:\n")
        for subject in subjects_PSIR:
            f.write("{}\n".format(subject))
        f.write("-----------------------------------\n")
        f.write("Subjects with STIR at every time point they have:\n")
        for subject in subjects_STIR:
            f.write("{}\n".format(subject))
        f.write("-----------------------------------\n")
        f.write("Subjects with PSIR and STIR at every time point they have:\n")
        for subject in subjects_PSIR_STIR:
            f.write("{}\n".format(subject))
        f.write("-----------------------------------\n")
        f.write("Subjects with no PSIR and no STIR at every time point they have:\n")
        for subject in subjects_no_PSIR_no_STIR:
            f.write("{}\n".format(subject))
        f.write("-----------------------------------\n")
        f.write("Subjects with no PSIR and no STIR at least once:\n")
        for subject in subjects_no_PSIR_no_STIR_once:
            f.write("{}\n".format(subject))
        f.write("-----------------------------------\n")
        f.write("Healthy control with PSIR at every time point they have:\n")
        for subject in subjects_hc_psir:
            f.write("{}\n".format(subject))
        f.write("-----------------------------------\n")
        f.write("Healthy control with STIR at every time point they have:\n")
        for subject in subjects_hc_stir:
            f.write("{}\n".format(subject))
        f.write("-----------------------------------\n")

    return None


if __name__ == '__main__':
    
    main()


