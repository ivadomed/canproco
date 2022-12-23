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

# Drop subjects, see: https://github.com/ivadomed/canproco/issues/13
subjects_to_exclude = ['sub-cal091', 'sub-cal155', 'sub-mon066', 'sub-mon033', 'sub-edm165', 'sub-mon006']


# T2w C2-C3 CSA values as in spine-generic SciData paper (doi: 10.1038/s41597-021-00941-8)
# Average +/− standard deviation (SD)
spine_generic_values = {
    "cal": (72.16, 3.06),       # ge
    "van": (72.76, 2.47),       # philips
    "mon": (72.76, 2.47),       # philips
    "edm": (75.44, 3.83),       # siemens
    "tor": (75.44, 3.83)        # siemens
}


# x position of individual sites
site_x_axis = {
    "cal": 0,
    "van": 1,
    "mon": 2,
    "edm": 3,
    "tor": 4
}


site_to_vendor = {
    "cal": "Calgary\nGE Discovery MR750",
    "van": "Vancouver\nPhilips Ingenia",
    "mon": "Montreal\n Philips Ingenia",
    "edm": "Edmonton\n Siemens Prisma",
    "tor": "Toronto\nSiemens Skyra",
}


def get_parser():
    parser = argparse.ArgumentParser(
        description="Generate figure for T2w C2-C3 CSA. The figure is saved to the same folder as the input .csv file."
    )
    parser.add_argument(
        '-i',
        required=True,
        metavar='<file_path>',
        help="input .csv file with CSA values")
    parser.add_argument(
        '-participants-file',
        required=True,
        metavar='<file_path>',
        help="participants.tsv file (includes pathology and phenotype columns)")

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


def add_spine_generic_values_per_vendor(ax, site, shift_i=0.15, shift_j=0.45):
    """
    Add mean and SD spine-generic values represented by rectangle and dashed line, respectively
    :param ax:
    :param site:
    :param shift_i: left shift from boxplots
    :param shift_j: right shift from boxplots
    :return:
    """
    mean = spine_generic_values[site][0]
    std = spine_generic_values[site][1]
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


def create_rainplot(metric_pd, fname_fig):
    """
    Create a rainplot (box + strip + violin)
    https://github.com/pog87/PtitPrince/blob/master/tutorial_python/raincloud_tutorial_python.ipynb
    :param metric_pd:
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
                      rain_alpha=1,      # individual points transparency - https://github.com/pog87/PtitPrince/blob/master/ptitprince/PtitPrince.py
                      alpha=.7           # violin plot transparency
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
        add_spine_generic_values_per_vendor(ax, site, shift_i=0.17, shift_j=0.34)

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


def create_violinplot(metric_pd, fname_fig):
    """
    Create a violinplot
    :param metric_pd:
    :param fname_fig:
    :return:
    """
    fig, ax = plt.subplots(figsize=(21, 7))
    ax = sns.violinplot(data=metric_pd,
                        x='site',
                        y='MEAN(area)',
                        hue='phenotype',
                        inner='point',      # include individual dots to the violin plots
                        cut=0               # limit the violin range within the range of the observed data
                        )
    # Change transparency
    for violin, alpha in zip(ax.collections[::2], len(ax.collections[::2]) * [0.8]):
        violin.set_alpha(alpha)
    plt.title('C2-C3 spinal cord cross-sectional area (CSA) from T2w', fontsize=FONTSIZE)
    ax.set_ylabel('CSA [$mm^2$]', fontsize=FONTSIZE)
    ax.set_xlabel('')
    # Move grid to background (i.e. behind other elements)
    ax.set_axisbelow(True)
    # Add horizontal grid lines
    ax.yaxis.grid(True)

    # Modify x-ticks labels
    ax.set_xticklabels(site_to_vendor.values())
    # Increase size of xticks and yticks
    plt.setp(ax.xaxis.get_majorticklabels(), fontsize=FONTSIZE)
    plt.setp(ax.yaxis.get_majorticklabels(), fontsize=FONTSIZE)

    # Increase legend title and text
    plt.setp(ax.get_legend().get_title(), fontsize=FONTSIZE)
    plt.setp(ax.get_legend().get_texts(), fontsize=FONTSIZE)

    # Add mean and SD spine-generic values
    for site in site_to_vendor.keys():
        add_spine_generic_values_per_vendor(ax, site)

    plt.tight_layout()

    # save figure
    plt.savefig(fname_fig, dpi=200)
    plt.close()
    print(f'Created: {fname_fig}.\n')


def compute_anova(metric_pd):
    """
    Compute ANOVA among phenotypes
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


def compute_kruskal(metric_pd):
    """
    Compute Kruskal-Wallis H-test among phenotypes
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
        # Compute Kruskal-Wallis H-test
        fvalue, pvalue = stats.kruskal(metric_pd_site[metric_pd_site['phenotype'] == 'RRMS']['MEAN(area)'],
                                       metric_pd_site[metric_pd_site['phenotype'] == 'PPMS']['MEAN(area)'],
                                       metric_pd_site[metric_pd_site['phenotype'] == 'RIS']['MEAN(area)'],
                                       metric_pd_site[metric_pd_site['phenotype'] == 'HC']['MEAN(area)'])
        print(f'Kruskal-Wallis H-test p-value: {pvalue}\n')
        # Post hoc Conover’s test
        # https://scikit-posthocs.readthedocs.io/en/latest/tutorial/#non-parametric-anova-with-post-hoc-tests
        posthoc = sp.posthoc_conover(metric_pd_site, val_col='MEAN(area)', group_col='phenotype', p_adjust='holm')
        print(f'{posthoc}\n')


def main():
    parser = get_parser()
    args = parser.parse_args()

    if os.path.isfile(args.i):
        # Read CSV file as pandas DF (read only columns with filename (contains sub_ID) and metric mean value)
        metric_pd = pd.read_csv(args.i, usecols=['Filename', 'MEAN(area)'])
    else:
        raise FileNotFoundError(f'{args.i} not found')

    if os.path.isfile(args.participants_file):
        # Read participants.tsv file (includes pathology and phenotype columns)
        participants_pd = pd.read_csv(args.participants_file, sep='\t')
    else:
        raise FileNotFoundError(f'{args.participants_file} not found')

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

    # Merge pathology and phenotype columns to the dataframe with CSA values
    metric_pd = pd.merge(metric_pd, participants_pd[['participant_id', 'pathology', 'phenotype']], how='left',
                         left_on='subject_id', right_on='participant_id')

    # Drop subjects, see: https://github.com/ivadomed/canproco/issues/13
    for subject in subjects_to_exclude:
        metric_pd.drop(metric_pd.loc[metric_pd['subject_id'] ==subject].index, inplace=True)
        print(f'Dropping {subject}')

    # Replace n/a in phenotype by HC to allow sorting in violinplot
    metric_pd['phenotype'].fillna(metric_pd['pathology'], inplace=True)

    # Create violin plot
    fname_fig = args.i.replace('.csv', '_violinplot.png')
    create_violinplot(metric_pd, fname_fig)

    # Create rain plot
    fname_fig = args.i.replace('.csv', '_rainplot.png')
    create_rainplot(metric_pd, fname_fig)

    # Compute ANOVA among phenotypes
    compute_anova(metric_pd)
    # Kruskal-Wallis H-test among phenotypes
    compute_kruskal(metric_pd)

    # Compute median, mean, std, cov per site and phenotype
    statistic = metric_pd.groupby(['site', 'phenotype']).agg([np.median, np.mean, np.std, stats.variation])
    print(f'\nDescriptive statistics:\n{statistic}')


if __name__ == "__main__":
    main()
