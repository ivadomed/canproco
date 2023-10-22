"""
Generate Lesion Frequency Maps (LFM) from the lesion masks of individual subjects in the PAM50 space

The script requires the SCT conda environment to be activated:
    source ${SCT_DIR}/python/etc/profile.d/conda.sh
    conda activate venv_sct

Author: Jan Valosek
Inspired by https://github.com/neuropoly/lesion-mapping/blob/master/spinalcord/3_generate_LFM.py from Charley Gros
"""

import os
import argparse
import pandas as pd
import numpy as np

from spinalcordtoolbox.image import Image, zeros_like

# Lateral and ventral corticospinal (CST) tracts
TRACTS_LST = ['PAM50_atlas_04.nii.gz', 'PAM50_atlas_05.nii.gz', 'PAM50_atlas_22.nii.gz', 'PAM50_atlas_23.nii.gz']

# Site to contrast
SITE_DCT = {
        'cal': 'STIR',
        'mon': 'PSIR',
        'tor': 'PSIR',
        'van': 'PSIR',
        'edm': 'PSIR'
        }


def get_parser():
    """
    parser function
    """

    parser = argparse.ArgumentParser(
        description='Generate Lesion Frequency Maps (LFM) from the lesion masks of individual subjects in the PAM50 '
                    'space',
        prog=os.path.basename(__file__).strip('.py')
    )
    parser.add_argument(
        '-ifolder',
        metavar="<file>",
        required=True,
        type=str,
        help='Path to the folder with processed MRI data. Example: /data_processed'
    )
    parser.add_argument(
        '-participants-tsv',
        metavar="<file>",
        required=True,
        type=str,
        help='Path to the participants.tsv file containing participant_id and phenotype_M0 columns. '
             'Example: canproco/participants.tsv'
    )
    parser.add_argument(
        '-ofolder',
        metavar="<folder>",
        required=True,
        type=str,
        help='Path to the output folder where LFM will be saved. Example: /results'
    )

    return parser


def clean_LFM(fname_out, fname_cord, fname_lvl):
    """
    Clean the LFM by removing voxels outside of the cord and above and below the vert levels C1 and C7
    """
    img, cord, lvl = Image(fname_out), Image(fname_cord), Image(fname_lvl)
    cord_data, lvl_data = cord.data, lvl.data
    del cord, lvl

    img.data[np.where(cord_data == 0)] = 0
    z_top = np.max(list(set(np.where(lvl_data == 1)[2]))) + 1
    z_bottom = np.min(list(set(np.where(lvl_data == 7)[2])))
    img.data[:, :, :z_bottom] = 0
    img.data[:, :, z_top:] = 0
    img.data[np.isnan(img.data)]=0.0
    img.data[img.data>1.0]=0.0

    img.save(fname_out)
    del img


def initialise_sumFile(fname_out, fname_standard):
    """
    Initialise the sum file by copying the PAM50 file
    """
    img_out = zeros_like(Image(fname_standard))
    img_out.save(fname_out)
    del img_out


def add_mask(fname_new, fname_out):
    """
    Add a mask to the sum file
    """
    img_new, img_in = Image(fname_new), Image(fname_out)
    img_out = zeros_like(img_in)
    img_out.data = img_new.data + img_in.data
    del img_new, img_in
    img_out.save(fname_out)
    del img_out


def mask_CST(fname_LFM, fname_LFM_CST, mask_lst):
    """
    Mask the LFM with the CST mask
    """
    img_lfm = Image(fname_LFM)
    img_cst = zeros_like(img_lfm)
    img_cst.data = img_lfm.data
    del img_lfm

    cst_mask_data = np.sum([Image(mask_fname).data for mask_fname in mask_lst], axis=0)
    cst_mask_data = (cst_mask_data > 0.0).astype(np.int_)

    img_cst.data[np.where(cst_mask_data == 0.0)] = 0.0
    img_cst.save(fname_LFM_CST)


def generate_LFM(df, fname_out, fname_out_cst, path_data):
    """
    Generate the LFM (Lesion Frequency Map)
    :param df: dataframe with participant_id and institution_id_M0 columns
    :param fname_out: output file name
    :param fname_out_cst: output file name for CST
    :param path_data: path to the data
    """
    path_pam50 = os.path.join(os.environ.get('SCT_DIR'), 'data/PAM50')
    pam50_cord = os.path.join(path_pam50, 'template', 'PAM50_cord.nii.gz')
    pam50_lvl = os.path.join(path_pam50, 'template', 'PAM50_levels.nii.gz')

    fname_out_lesion = fname_out.split('_LFM.nii.gz')[0] + '_sumLesion.nii.gz'
    fname_out_cord = fname_out.split('_LFM.nii.gz')[0] + '_sumCord.nii.gz'
    initialise_sumFile(fname_out_lesion, pam50_cord)
    initialise_sumFile(fname_out_cord, pam50_cord)

    for index, row in df.iterrows():
        participant_id = row.participant_id

        path_subject = os.path.join(path_data, participant_id, 'ses-M0', 'anat')
        # Lesion in PAM50 space
        lesion_path = os.path.join(path_subject,
                                   f'{participant_id}_ses-M0_{SITE_DCT[row.institution_id_M0]}_lesion-manual_bin_reg.nii.gz')
        # Spinal cord in PAM50 space
        cord_path = os.path.join(path_subject,
                                 f'{participant_id}_ses-M0_{SITE_DCT[row.institution_id_M0]}_seg-manual_reg.nii.gz')

        if os.path.isfile(lesion_path) and os.path.isfile(cord_path):
            print(participant_id)
            add_mask(lesion_path, fname_out_lesion)
            add_mask(cord_path, fname_out_cord)

    os.system('sct_maths'
              ' -i ' + fname_out_lesion +
              ' -div ' + fname_out_cord +
              ' -o ' + fname_out)

    clean_LFM(fname_out, pam50_cord, pam50_lvl)
    mask_CST(fname_out, fname_out_cst, [os.path.join(path_pam50, 'atlas', t) for t in TRACTS_LST])


def main():
    # Parse the command line arguments
    parser = get_parser()
    args = parser.parse_args()

    # Path to the folder containing the processed data
    path_folder = args.ifolder
    # Path to the output folder
    path_out = args.ofolder
    if not os.path.isdir(path_out):
        os.makedirs(path_out)

    path_participants_tsv = args.participants_tsv
    # Read the participants.tsv file
    participants_df = pd.read_csv(path_participants_tsv, sep='\t',
                                  usecols=['participant_id', 'institution_id_M0', 'phenotype_M0', 'edss_M0'])

    # Loop across the subgroups
    for subgroup in ['all', 'RRMS', 'PPMS', 'RIS', 'edss_low', 'edss_med', 'edss_high']:
        path_lfm = os.path.join(path_out, 'spinalcord_LFM_' + subgroup + '.nii.gz')
        path_lfm_cst = os.path.join(path_out, 'spinalcord_LFM_CST_' + subgroup + '.nii.gz')
        if subgroup == 'all':
            lfm_df = participants_df
        elif subgroup in ['RRMS', 'PPMS', 'RIS']:
            lfm_df = participants_df[participants_df.phenotype_M0 == subgroup]
        elif subgroup.startswith('edss_'):
            if subgroup.endswith('low'):
                lfm_df = participants_df[participants_df.edss_M0 <= 2.5]
            elif subgroup.endswith('high'):
                lfm_df = participants_df[participants_df.edss_M0 >= 6.0]
            elif subgroup.endswith('med'):
                lfm_df = participants_df[(participants_df.edss_M0 < 6.0) & (participants_df.edss_M0 > 2.5)]

        if not os.path.isfile(path_lfm) or not os.path.isfile(path_lfm_cst):
            print(f'\nGenerating the LFM with {subgroup} subjects ({str(len(lfm_df.index))}).')
            generate_LFM(lfm_df, path_lfm, path_lfm_cst, path_folder)


if __name__ == "__main__":
    main()
