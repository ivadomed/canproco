#!/usr/bin/env python
#
# Generate figure for T2w C2-C3 CSA
#
# Authors: Jan Valosek, Julien Cohen-Adad

import os
import re
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import ptitprince as pt

from scipy import stats
import scikit_posthocs as sp
import matplotlib.patches as patches
from matplotlib.cm import get_cmap
from matplotlib.lines import Line2D

FONTSIZE=18

# Drop canproco subjects, see: https://github.com/ivadomed/canproco/issues/13
subjects_to_exclude_canproco = ['sub-cal091', 'sub-cal155', 'sub-mon066', 'sub-mon033', 'sub-edm165', 'sub-mon006']

# Drop spine-generic subjects
# site beijingVerio - different TR and FA causing biases in the segmentation volume; see doi:10.1038/s41597-021-00941-8
subjects_to_exclude_spinegeneric = ['sub-beijingVerio01', 'sub-beijingVerio02', 'sub-beijingVerio03',
                                    'sub-beijingVerio04']

# x position of individual sites
site_x_axis = {
    "cal": 0,
    "van": 1,
    "mon": 2,
    "edm": 3,
    "tor": 4
}

# xticks labels
site_to_vendor = {
    "cal": "Calgary\nGE Discovery MR750",
    "van": "Vancouver\nPhilips Ingenia",
    "mon": "Montreal\n Philips Ingenia",
    "edm": "Edmonton\n Siemens Prisma",
    "tor": "Toronto\nSiemens Skyra",
}

site_to_manufacturer = {
    "cal": "GE",
    "van": "Philips",
    "mon": "Philips",
    "edm": "Siemens",
    "tor": "Siemens",
}


def get_parser():
    parser = argparse.ArgumentParser(
        description="Generate figure for T2w C2-C3 CSA. The figure is saved to the same folder as the input .csv file."
    )
    parser.add_argument(
        '-i-canproco',
        required=True,
        metavar='<file_path>',
        help="input .csv file with canproco CSA values")
    parser.add_argument(
        '-i-spinegeneric',
        required=True,
        metavar='<file_path>',
        help="input .csv file with spine-generic CSA values")
    parser.add_argument(
        '-participants-file-canproco',
        required=True,
        metavar='<file_path>',
        help="canproco participants.tsv file (includes pathology and phenotype columns)")
    parser.add_argument(
        '-participants-file-spinegeneric',
        required=True,
        metavar='<file_path>',
        help="spine-generic participants.tsv file (includes vendor column)")

    return parser


# TODO - this function could be merged with function in utils.py
def fetch_subject_and_site(filename_path):
    """
    Get subject ID and site from the input BIDS-compatible filename or file path
    The function works both on absolute file path as well as filename
    :param filename_path: input nifti filename (e.g., sub-001_ses-01_T1w.nii.gz) or file path
    (e.g., /home/user/MRI/bids/derivatives/labels/sub-001/ses-01/anat/sub-001_ses-01_T1w.nii.gz
    :return: subject_id: subject ID (e.g., sub-001)
    :return: sessionID: session ID (e.g., ses-01)
    :return: filename: nii filename (e.g., sub-001_ses-01_T1w.nii.gz)
    """

    _, filename = os.path.split(filename_path)              # Get just the filename (i.e., remove the path)
    subject = re.search('sub-(.*?)[_/]', filename_path)
    subject_id = subject.group(0)[:-1] if subject else ""    # [:-1] removes the last underscore or slash
    site = subject_id.split('-')[1][:3]
    # REGEX explanation
    # \d - digit
    # \d? - no or one occurrence of digit
    # *? - match the previous element as few times as possible (zero or more times)

    return subject_id, site


def add_spine_generic_values_per_vendor(ax, site, spinegeneric_pd, shift_i=0.15, shift_j=0.45):
    """
    Add mean and SD spine-generic values represented by rectangle and dashed line, respectively
    :param ax:
    :param site:
    :param spinegeneric_pd: Pandas DataFrame with spine-generic CSA values
    :param shift_i: left shift from boxplots
    :param shift_j: right shift from boxplots
    :return:
    """
    # compute mean within vendor (mean of the within-site means)
    # https://github.com/spine-generic/spine-generic/blob/master/spinegeneric/cli/generate_figure.py#L423
    mean = float(spinegeneric_pd[spinegeneric_pd['manufacturer'] == site_to_manufacturer[site]].groupby(['site']).agg([np.mean]).mean())
    # compute std within vendor (std of the within-site means)
    # https://github.com/spine-generic/spine-generic/blob/master/spinegeneric/cli/generate_figure.py#L425
    std = float(spinegeneric_pd[spinegeneric_pd['manufacturer'] == site_to_manufacturer[site]].groupby(['site']).agg([np.mean]).std())
    print(f'{site_to_manufacturer[site]}: {mean} +- {std}')
    x = site_x_axis[site]
    # add rectangle for variance
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Rectangle.html#matplotlib-patches-rectangle
    rect = patches.Rectangle(xy=(x - shift_i, mean - std),      # shift_i is used to cover boxplots (but not violinpltos)
                             width=shift_j,                     # shift_j is used to cover also individual points (when 'move' rainplot param is used)
                             height=2*std,
                             edgecolor=None,
                             facecolor='gray',
                             alpha=0.3)
    ax.add_patch(rect)
    # add dashed line for mean value
    ax.plot([x - shift_i, x + shift_j-shift_i], [mean, mean], "k--", alpha=0.3)


def create_rainplot(metric_pd, spinegeneric_pd, fname_fig):
    """
    Create a rainplot (box + strip + violin)
    https://github.com/pog87/PtitPrince/blob/master/tutorial_python/raincloud_tutorial_python.ipynb
    :param metric_pd: Pandas DataFrame with canproco CSA values
    :param spinegeneric_pd: Pandas DataFrame with spine-generic CSA values
    :param fname_fig:
    :return:
    """
    fig, ax = plt.subplots(figsize=(21, 7))
    ax = pt.RainCloud(data=metric_pd,
                      x='site',
                      y='MEAN(area)',
                      hue='phenotype',
                      order=site_to_vendor.keys(),
                      palette="Set1",
                      linewidth=0,       # violionplot border line (0 - no line)
                      width_viol=.5,     # violinplot width
                      width_box=.3,      # boxplot width
                      dodge=True,        # move boxplots next to each other
                      move=0,            # move individual observations next to the boxplots (0 - no move)
                      rain_alpha=1,      # individual points transparency - https://github.com/pog87/PtitPrince/blob/23debd9b70fca94724a06e72e049721426235f50/ptitprince/PtitPrince.py#L707
                      alpha=.7,          # violin plot transparency
                      box_showmeans=True,   # show mean value inside the boxplots
                      box_meanprops={'marker': '^', 'markerfacecolor': 'black', 'markeredgecolor': 'black',
                                     'markersize': '6'}
                      )
    plt.title('C2-C3 spinal cord cross-sectional area (CSA) from T2w', fontsize=FONTSIZE)
    ax.set_ylabel('CSA [$mm^2$]', fontsize=FONTSIZE)
    ax.set_xlabel('')
    # Move grid to background (i.e. behind other elements)
    ax.set_axisbelow(True)
    # Add horizontal grid lines
    ax.yaxis.grid(True)

    # Modify x-ticks labels (set full site name and scanner type)
    ax.set_xticklabels(site_to_vendor.values())
    # Increase size of xticks and yticks
    plt.setp(ax.xaxis.get_majorticklabels(), fontsize=FONTSIZE)
    plt.setp(ax.yaxis.get_majorticklabels(), fontsize=FONTSIZE)

    # Increase legend title and text
    plt.setp(ax.get_legend().get_title(), fontsize=FONTSIZE)
    plt.setp(ax.get_legend().get_texts(), fontsize=FONTSIZE)

    # Change boxplot opacity (.0 means transparent)
    # https://github.com/mwaskom/seaborn/issues/979#issuecomment-1144615001
    for patch in ax.patches:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .0))

    # Add mean and SD spine-generic values
    for site in site_to_vendor.keys():
        add_spine_generic_values_per_vendor(ax, site, spinegeneric_pd,  shift_i=0.17, shift_j=0.34)

    # LEGEND
    _, labels = ax.get_legend_handles_labels()
    n_plots = 3     # violinplot + boxplot + pointplot = 3
    # create colorful lines
    lines = [Line2D([0], [0], color=value, linestyle='-', linewidth=10)
             for value in get_cmap('Set1').colors[:4]]
    # crete a line for spine-generic
    line_spine_generic = Line2D([0], [0], color='gray', linestyle='--', linewidth=3)
    lines.append(line_spine_generic)
    # legend labels
    legend_labels = labels[0:len(labels) // n_plots]
    legend_labels.append(u'spine-generic\nmean \u00B1 SD')
    # Move legend closer to the plot (bbox_to_anchor) and set the length of the '-' (handlelength)
    legend = plt.legend(lines, legend_labels, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.,
                        handlelength=2, title='Phenotype', fontsize=FONTSIZE-3, title_fontsize=FONTSIZE)
    # Change box's frame color to black
    frame = legend.get_frame()
    frame.set_edgecolor('black')

    plt.tight_layout()

    # save figure
    plt.savefig(fname_fig, dpi=200)
    plt.close()
    print(f'Created: {fname_fig}.\n')


def compute_anova_per_site(metric_pd):
    """
    Compute ANOVA among phenotypes persite
    :param metric_pd:
    :return:
    """
    # Loop across sites
    for site in site_to_vendor.keys():
        # Get values only for given site
        metric_pd_site = metric_pd[metric_pd['site'] == site]
        # Compute one-way ANOVA
        print(f'{site}')
        print(metric_pd_site.groupby(['phenotype']).size())
        fvalue, pvalue = stats.f_oneway(metric_pd_site[metric_pd_site['phenotype'] == 'RRMS']['MEAN(area)'],
                                        metric_pd_site[metric_pd_site['phenotype'] == 'PPMS']['MEAN(area)'],
                                        metric_pd_site[metric_pd_site['phenotype'] == 'RIS']['MEAN(area)'],
                                        metric_pd_site[metric_pd_site['phenotype'] == 'HC']['MEAN(area)'])
        print(f'ANOVA p-value: {pvalue}')


def compute_kruskal_per_site(metric_pd):
    """
    Compute Kruskal-Wallis H-test among phenotypes persite
    :param metric_pd:
    :return:
    """
    # Loop across sites
    for site in site_to_vendor.keys():
        # Get values only for given site
        metric_pd_site = metric_pd[metric_pd['site'] == site]
        # Compute Kruskal-Wallis H-test
        print(f'{site}')
        print(metric_pd_site.groupby(['phenotype']).size())
        fvalue, pvalue = stats.kruskal(metric_pd_site[metric_pd_site['phenotype'] == 'RRMS']['MEAN(area)'],
                                       metric_pd_site[metric_pd_site['phenotype'] == 'PPMS']['MEAN(area)'],
                                       metric_pd_site[metric_pd_site['phenotype'] == 'RIS']['MEAN(area)'],
                                       metric_pd_site[metric_pd_site['phenotype'] == 'HC']['MEAN(area)'])
        print(f'Kruskal-Wallis H-test p-value: {pvalue}\n')
        # Post hoc Conover’s test
        # https://scikit-posthocs.readthedocs.io/en/latest/tutorial/#non-parametric-anova-with-post-hoc-tests
        posthoc = sp.posthoc_conover(metric_pd_site, val_col='MEAN(area)', group_col='phenotype', p_adjust='holm')
        print(f'{posthoc}\n')


def compute_kruskal(metric_pd):
    """
    Compute Kruskal-Wallis H-test among phenotypes on the whole cohort
    :param metric_pd:
    :return:
    """
    # Compute Kruskal-Wallis H-test
    print(f'Whole cohort')
    print(metric_pd.groupby(['phenotype']).size())
    fvalue, pvalue = stats.kruskal(metric_pd[metric_pd['phenotype'] == 'RRMS']['MEAN(area)'],
                                   metric_pd[metric_pd['phenotype'] == 'PPMS']['MEAN(area)'],
                                   metric_pd[metric_pd['phenotype'] == 'RIS']['MEAN(area)'],
                                   metric_pd[metric_pd['phenotype'] == 'HC']['MEAN(area)'])
    print(f'Kruskal-Wallis H-test p-value: {pvalue}\n')
    # Post hoc Conover’s test
    # https://scikit-posthocs.readthedocs.io/en/latest/tutorial/#non-parametric-anova-with-post-hoc-tests
    posthoc = sp.posthoc_conover(metric_pd, val_col='MEAN(area)', group_col='phenotype', p_adjust='holm')
    print(f'{posthoc}\n')


def compare_healthy_controls(canproco_pd, spinegeneric_pd):
    """
    Compare CSA values between canproco healthy controls and spine-generic per manufacturer. For example, Calgary site
    (GE) is compared to spine-generic GE values.
    :param canproco_pd:
    :param spinegeneric_pd:
    :return:
    """
    # Loop across sites
    for site, manufacturer in site_to_manufacturer.items():
        # Get canproco values for given site
        canproco_values = canproco_pd[canproco_pd['site'] == site]['MEAN(area)']
        # Get spine-generic values for given manufacturer
        spinegeneric_values = spinegeneric_pd[spinegeneric_pd['manufacturer'] == manufacturer]['MEAN(area)']
        # Perform the Mann-Whitney U rank test
        _, pvalue = stats.mannwhitneyu(canproco_values, spinegeneric_values)
        print(f'canproco site {site}, spine-generic manufacturer {manufacturer}: {pvalue}')


def read_csv_file(file_path, subjects_to_exclude, columns_to_read=['Filename', 'MEAN(area)']):
    """
    Read .csv file generated by the sct_process_segmentation. Then, fetch subjectID and site and them as new columns.
    :param file_path:
    :param subjects_to_exclude: list of subjects to exlude
    :param columns_to_read: columns to read
    :return: Pandas DataFrame
    """
    # Check if file exists
    if os.path.isfile(file_path):
        # Read CSV file as pandas DF (read only columns with filename (contains sub_ID) and metric mean value)
        metric_pd = pd.read_csv(file_path, usecols=columns_to_read)
    else:
        raise FileNotFoundError(f'{file_path} not found')

    # Fetch subjectID and site for each subject and insert it to the metric_pd
    site_list = list()
    subject_id_list = list()
    # Loop across rows
    for index, row in metric_pd.iterrows():
        # Fetch subjectID and site
        subject_id, site = fetch_subject_and_site(row['Filename'])
        subject_id_list.append(subject_id)
        site_list.append(site)
    # Insert subjectID and site lists into pandas DF
    metric_pd.insert(1, "subject_id", subject_id_list)
    metric_pd.insert(2, "site", site_list)

    # Drop subjects, see: https://github.com/ivadomed/canproco/issues/13
    for subject in subjects_to_exclude:
        metric_pd.drop(metric_pd.loc[metric_pd['subject_id'] == subject].index, inplace=True)
        print(f'Dropping {subject}')

    return metric_pd


def read_participants_file(file_path):
    """
    Read participants.tsv file and return Pandas DataFrame
    :param file_path:
    :return:
    """
    if os.path.isfile(file_path):
        participants_pd = pd.read_csv(file_path, sep='\t')
        return participants_pd
    else:
        raise FileNotFoundError(f'{file_path} not found')


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Read .csv file for canproco subjects
    canproco_pd = read_csv_file(args.i_canproco, subjects_to_exclude_canproco)
    # Read canproco participants.tsv file (includes pathology and phenotype columns)
    canproco_participants_pd = read_participants_file(args.participants_file_canproco)

    # Read .csv file for spine-generic subjects
    spinegeneric_pd = read_csv_file(args.i_spinegeneric, subjects_to_exclude_spinegeneric)
    # Read spine-generic participants.tsv file (includes manufacturer column)
    spinegeneric_participants_pd = read_participants_file(args.participants_file_spinegeneric)

    # Merge pathology and phenotype columns to the canproco dataframe with CSA values
    canproco_pd = pd.merge(canproco_pd, canproco_participants_pd[['participant_id', 'pathology', 'phenotype']],
                           how='left', left_on='subject_id', right_on='participant_id')

    # Replace n/a in phenotype by HC to allow sorting in violinplot
    canproco_pd['phenotype'].fillna(canproco_pd['pathology'], inplace=True)

    # Merge manufacturer column to the spine-generic dataframe with CSA values
    spinegeneric_pd = pd.merge(spinegeneric_pd, spinegeneric_participants_pd[['participant_id', 'manufacturer']],
                               how='left', left_on='subject_id', right_on='participant_id')

    # Create rain plot
    fname_fig = args.i_canproco.replace('.csv', '_rainplot.png')
    create_rainplot(canproco_pd, spinegeneric_pd, fname_fig)

    # Compute ANOVA among phenotypes
    compute_anova_per_site(canproco_pd)
    # Kruskal-Wallis H-test among phenotypes
    compute_kruskal_per_site(canproco_pd)
    # Compute Kruskal-Wallis H-test among phenotypes on the whole cohort
    compute_kruskal(canproco_pd)

    # Compute median, mean, std, cov persite and phenotype
    statistic = canproco_pd.groupby(['site', 'phenotype']).agg([np.median, np.mean, np.std, stats.variation])
    print(f'\nDescriptive statistics:\n{statistic}')

    # Compute median, mean, std, cov per phenotype on the whole cohort
    statistic = canproco_pd.groupby(['phenotype']).agg([np.median, np.mean, np.std, stats.variation])
    print(f'\nDescriptive statistics:\n{statistic}')

    # Compare CSA values between canproco healthy controls and spine-generic per manufacturer
    compare_healthy_controls(canproco_pd, spinegeneric_pd)


if __name__ == "__main__":
    main()
