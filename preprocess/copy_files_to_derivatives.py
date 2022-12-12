#!/usr/bin/env python
#
# Copy manually corrected files (segmentations, vertebral labeling, etc.) from the source preprocessed dataset to the
# git-annex BIDS dataset's derivatives folder
#
# Authors: Jan Valosek

import argparse
import glob
import os
import shutil
import utils
import re


def get_parser():
    """
    parser function
    """

    parser = argparse.ArgumentParser(
        description='Copy manually corrected files (segmentations, vertebral labeling, etc.) from the source '
                    'preprocessed dataset to the git-annex BIDS derivatives folder',
        formatter_class=utils.SmartFormatter,
        prog=os.path.basename(__file__).strip('.py')
    )
    parser.add_argument(
        '-path-in',
        metavar="<folder>",
        required=True,
        type=str,
        help='Path to the folder with manually corrected files (usually derivatives). The script assumes that labels '
             'folder is located in the provided folder.'
    )
    parser.add_argument(
        '-path-out',
        metavar="<folder>",
        required=True,
        type=str,
        help='Path to the BIDS dataset where manually corrected files will be copied. Files will be copied to the '
             'derivatives/label folder.'
    )

    return parser


# TODO - merge this function with function in utils.py
def fetch_subject_and_session(filename_path):
    """
    Get subject ID, session ID and filename from the input BIDS-compatible filename or file path
    The function works both on absolute file path as well as filename
    :param filename_path: input nifti filename (e.g., sub-001_ses-01_T1w.nii.gz) or file path
    (e.g., /home/user/MRI/bids/derivatives/labels/sub-001/ses-01/anat/sub-001_ses-01_T1w.nii.gz
    :return: subjectID: subject ID (e.g., sub-001)
    :return: sessionID: session ID (e.g., ses-01)
    :return: filename: nii filename (e.g., sub-001_ses-01_T1w.nii.gz)
    """

    _, filename = os.path.split(filename_path)              # Get just the filename (i.e., remove the path)
    subject = re.search('sub-(.*?)[_/]', filename_path)
    subjectID = subject.group(0)[:-1] if subject else ""    # [:-1] removes the last underscore or slash
    session = re.findall(r'ses-..', filename_path)
    sessionID = session[0] if session else ""               # Return None if there is no session
    # REGEX explanation
    # \d - digit
    # \d? - no or one occurrence of digit
    # *? - match the previous element as few times as possible (zero or more times)

    return subjectID, sessionID, filename


def main():

    # Parse the command line arguments
    parser = get_parser()
    args = parser.parse_args()

    # Check if path_in exists
    if os.path.isdir(args.path_in):
        path_in = os.path.join(os.path.abspath(args.path_in), 'labels')
    else:
        raise NotADirectoryError(f'{args.path_in} does not exist.')

    # Check if path_out exists
    if os.path.isdir(args.path_out):
        path_out = os.path.join(os.path.abspath(args.path_out), 'labels')
    else:
        raise NotADirectoryError(f'{args.path_out} does not exist.')

    # Loop across files in input dataset
    for path_file_in in sorted(glob.glob(path_in + '/**/*.nii.gz', recursive=True)):
        sub, ses, filename = fetch_subject_and_session(path_file_in)
        # Construct path for the output file
        path_file_out = os.path.join(path_out, sub, ses, filename)
        # Check if subject's folder exists in the output dataset, if not, create it
        path_subject_folder_out = os.path.join(path_out, sub, ses)
        if not os.path.isdir(path_subject_folder_out):
            os.makedirs(path_subject_folder_out)
            print(f'Creating directory: {path_subject_folder_out}')
        # Copy nii and json files to the output dataset
        # TODO - consider rsync instead of shutil.copy
        shutil.copy(path_file_in, path_file_out)
        print(f'Copying: {path_file_in} to {path_file_out}')
        path_file_json_in = path_file_in.replace('nii.gz', 'json')
        path_file_json_out = path_file_out.replace('nii.gz', 'json')
        shutil.copy(path_file_json_in, path_file_json_out)
        print(f'Copying: {path_file_json_in} to {path_file_json_out}')


if __name__ == '__main__':
    main()
