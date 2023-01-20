#!/usr/bin/env python
#
# The script performs:
#   - generate rainplots for T2w C2-C3 canproco CSA per site and includes also spine-generic values
#   - compute Kruskal-Wallis H-test among phenotypes persite and on the whole cohort
#   - compare CSA values between canproco healthy controls and spine-generic per manufacturer
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
import pingouin as pg

from scipy import stats
import scikit_posthocs as sp
import matplotlib.patches as patches
from matplotlib.cm import get_cmap
from matplotlib.lines import Line2D
from sklearn.linear_model import LinearRegression

FONTSIZE=18
FONTSIZE_CORR=25

colormap = 'Set1'
color_pallete = get_cmap(colormap).colors[:4]

# Drop canproco subjects, see: https://github.com/ivadomed/canproco/issues/13
# 'sub-van175' - pending EDSS
subjects_to_exclude_canproco = ['sub-cal091', 'sub-cal155', 'sub-mon066', 'sub-mon033', 'sub-edm165', 'sub-mon006',
                                'sub-van175']

# Drop spine-generic subjects
# site beijingVerio - different TR and FA causing biases in the segmentation volume; see doi:10.1038/s41597-021-00941-8
subjects_to_exclude_spinegeneric = ['sub-beijingVerio01', 'sub-beijingVerio02', 'sub-beijingVerio03',
                                    'sub-beijingVerio04']

# x position of individual sites
site_x_axis = {
    "all": 0,
    "cal": 1,
    "van": 2,
    "mon": 3,
    "edm": 4,
    "tor": 5
}

# xticks labels
site_to_vendor = {
    "all": "All sites",
    "cal": "Calgary\nGE Discovery MR750",
    "van": "Vancouver\nPhilips Ingenia",
    "mon": "Montreal\nPhilips Ingenia",
    "edm": "Edmonton\nSiemens Prisma",
    "tor": "Toronto\nSiemens Skyra",
}

# title for EDSS vs CSA correlation figure
site_to_vendor_title = {
    "all": "All sites",
    "cal": "Calgary - GE Discovery MR750",
    "van": "Vancouver - Philips Ingenia",
    "mon": "Montreal - Philips Ingenia",
    "edm": "Edmonton - Siemens Prisma",
    "tor": "Toronto - Siemens Skyra",
}

site_to_manufacturer = {
    "all": "All sites",
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
    Add mean and SD spine-generic values represented by rectangle and dashed line, respectively.
    Mean and SD values are added per-vendor (Philips, GE, or Siemens), depending on the MR vendor for given site.
    :param ax:
    :param site: currently processed site
    :param spinegeneric_pd: Pandas DataFrame with spine-generic CSA values
    :param shift_i: left shift from boxplots
    :param shift_j: right shift from boxplots
    :return:
    """
    if site == 'all':
        # compute mean within vendor (mean of the within-site means)
        # https://github.com/spine-generic/spine-generic/blob/master/spinegeneric/cli/generate_figure.py#L423
        mean = float(spinegeneric_pd.groupby(['site']).agg([np.mean]).mean())
        # compute std within vendor (std of the within-site means)
        # https://github.com/spine-generic/spine-generic/blob/master/spinegeneric/cli/generate_figure.py#L425
        std = float(spinegeneric_pd.groupby(['site']).agg([np.mean]).std())
    else:
        # compute mean within vendor (mean of the within-site means)
        # https://github.com/spine-generic/spine-generic/blob/master/spinegeneric/cli/generate_figure.py#L423
        mean = float(spinegeneric_pd[spinegeneric_pd['manufacturer'] == site_to_manufacturer[site]].groupby(['site']).
                     agg([np.mean]).mean())
        # compute std within vendor (std of the within-site means)
        # https://github.com/spine-generic/spine-generic/blob/master/spinegeneric/cli/generate_figure.py#L425
        std = float(spinegeneric_pd[spinegeneric_pd['manufacturer'] == site_to_manufacturer[site]].groupby(['site']).
                    agg([np.mean]).std())

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
    # duplicated the whole dataframe, but change site to 'all' --> this will allow to add 'All sites' to rainplot
    temp_pd = metric_pd.copy()
    temp_pd['site'] = 'all'
    concat_pd = pd.concat([metric_pd, temp_pd])

    fig, ax = plt.subplots(figsize=(21, 7))
    ax = pt.RainCloud(data=concat_pd,
                      x='site',
                      y='MEAN(area)',
                      hue='phenotype',
                      order=site_to_vendor.keys(),
                      palette=colormap,
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
    # Set ylim to have space for markers of significance
    ax.set_ylim([50, 105])

    ax.set_xlabel('')
    # Move grid to background (i.e. behind other elements)
    ax.set_axisbelow(True)
    # Add horizontal grid lines
    ax.yaxis.grid(True)

    # Modify x-ticks labels (set full site name and scanner type)
    ax.set_xticklabels(site_to_vendor.values())
    # Make 'All sites' label bold
    label_all_sites = ax.get_xticklabels()[0]
    label_all_sites.set_fontweight('bold')

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

    #plt.text(0, 50, '20', fontsize=10, horizontalalignment='center')

    # Add mean and SD spine-generic values
    for site in site_to_vendor.keys():
        add_spine_generic_values_per_vendor(ax, site, spinegeneric_pd,  shift_i=0.17, shift_j=0.34)

    # LEGEND
    _, labels = ax.get_legend_handles_labels()
    n_plots = 3     # violinplot + boxplot + pointplot = 3
    # create colorful lines
    lines = [Line2D([0], [0], color=value, linestyle='-', linewidth=10)
             for value in color_pallete]
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


def format_pvalue(p_value, alpha=0.001, decimal_places=3, include_space=False, include_equal=True):
    """
    Format p-value.
    If the p-value is lower than alpha, format it to "<0.001", otherwise, round it to three decimals
    :param p_value: input p-value as a float
    :param alpha: significance level
    :param decimal_places: number of decimal places the p-value will be rounded
    :param include_space: include space or not (e.g., ' = 0.06')
    :param include_equal: include equal sign ('=') to the p-value (e.g., '=0.06') or not (e.g., '0.06')
    :return: p_value: the formatted p-value (e.g., '<0.05') as a str
    """
    if include_space:
        space = ' '
    else:
        space = ''

    # If the p-value is lower than alpha, return '<alpha' (e.g., <0.001)
    if p_value < alpha:
        p_value = space + "<" + space + str(alpha)
    # If the p-value is greater than alpha, round it number of decimals specified by decimal_places
    else:
        if include_equal:
            p_value = space + '=' + space + str(round(p_value, decimal_places))
        else:
            p_value = space + str(round(p_value, decimal_places))

    return p_value


def compute_partial_correlation(canproco_pd, site):
    """
    Compute partial correlation with phenotype as a covariate
    :param canproco_pd:
    :param site:
    :return:
    """
    # Work only with MS patients
    if site == 'all':
        ms_pd = canproco_pd[canproco_pd['pathology'] == 'MS']
    else:
        ms_pd = canproco_pd[(canproco_pd['pathology'] == 'MS') & (canproco_pd['site'] == site)]
    # Convert str to int (to be compatible with partial correlation)
    ms_pd = ms_pd.replace({'phenotype': {'RRMS': 0, 'PPMS': 1, 'RIS': 2}})
    stats = pg.partial_corr(data=ms_pd, x='MEAN(area)', y='edss_M0', covar='phenotype', method='spearman')
    r = float(stats['r'])
    p_val = float(stats['p-val'])

    return r, p_val


def compute_correlation(csa, edss):

    r, p_val = stats.pearsonr(csa, edss)

    return r, p_val


def compute_regression(x, y):
    """
    Compute a linear regression between x and y:
    y = Slope * x + Intercept
    https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LinearRegression.html
    :param x: ndarray: input - regressor
    :param y: ndarray: output - response
    :return x_vals: ndarray: x values for the linear fit plot
    :return y_vals: ndarray: y values for the linear fit plot
    """
    # Make sure we are working with numpy arrays
    if isinstance(x, pd.Series):
        x = x.to_numpy()
    if isinstance(y, pd.Series):
        y = y.to_numpy()

    # Create an instance of the class LinearRegression, which will represent the regression model
    linear_regression = LinearRegression()
    # Perform linear regression (compute slope and intercept)
    linear_regression.fit(x.reshape(-1, 1), y.reshape(-1, 1))
    intercept = linear_regression.intercept_        # underscore indicates that an attribute is estimated
    slope = linear_regression.coef_                 # underscore indicates that an attribute is estimated

    # Get x and y values to plot the linear fit
    x_vals = np.array([x.min(), x.max()])
    y_vals = intercept + slope * x_vals
    y_vals = np.squeeze(y_vals)                     # change shape from (1,N) to (N,)

    return x_vals, y_vals


def create_csa_edss_correlation_figure_persite(canproco_pd, fname_fig):
    """
    Plot the relationship between EDSS score and CSA per-site and per-phenotype. Also, plot linear fit per-phenotype and
    for the whole cohort.
    :param canproco_pd:
    :param fname_fig:
    :return:
    """
    # Create main figure
    fig, axes = plt.subplots(2, 3, figsize=(20, 14), sharey=True)
    # Flatten 2D array into 1D to allow iteration by loop
    ax = axes.ravel()
    # Loop across sites (all means all sites together)
    for index, site in enumerate(site_to_vendor_title.keys()):
        # Compute partial correlation (with phenotype as a covariate)
        r, p_val = compute_partial_correlation(canproco_pd, site)
        print(f'{site}: Partial correlation EDSS vs CSA: r={r}, p-value{format_pvalue(p_val, alpha=0.05)}')
        # Compute linear regression for all MS patients together (i.e., across all phenotypes) --> ['pathology'] == 'MS'
        if site == 'all':
            csa = canproco_pd[canproco_pd['pathology'] == 'MS']['MEAN(area)']
            edss = canproco_pd[canproco_pd['pathology'] == 'MS']['edss_M0']
            #phen = canproco_pd[canproco_pd['pathology'] == 'MS']['phenotype']
        else:
            csa = canproco_pd[(canproco_pd['pathology'] == 'MS') & (canproco_pd['site'] == site)]['MEAN(area)']
            edss = canproco_pd[(canproco_pd['pathology'] == 'MS') & (canproco_pd['site'] == site)]['edss_M0']
            #phen = canproco_pd[(canproco_pd['pathology'] == 'MS') & (canproco_pd['site'] == site)]['phenotype']
        x_vals, y_vals = compute_regression(csa, edss)
        ax[index].plot(x_vals, y_vals, '--', color='black', alpha=.5, linewidth=3)

        # Insert text with corr coef and pval into every subplot/axis
        ax[index].annotate('r={}; p{}'.format(round(r, 2), format_pvalue(p_val, alpha=0.05)), xy=(.98, .9),
                           xycoords='axes fraction', fontsize=FONTSIZE_CORR-5, xytext=(-5, 5), textcoords='offset points',
                           ha='right', va='bottom', bbox=dict(edgecolor='black', facecolor='none', boxstyle='round'))

        for color, phenotype in enumerate(['RRMS', 'PPMS', 'RIS']):
            # Prepare variables for plotting
            if site == 'all':
                csa = canproco_pd[canproco_pd['phenotype'] == phenotype]['MEAN(area)']
                edss = canproco_pd[canproco_pd['phenotype'] == phenotype]['edss_M0']
                r, p_val = compute_correlation(csa, edss)
            else:
                csa = canproco_pd[(canproco_pd['phenotype'] == phenotype) & (canproco_pd['site'] == site)]['MEAN(area)']
                edss = canproco_pd[(canproco_pd['phenotype'] == phenotype) & (canproco_pd['site'] == site)]['edss_M0']
                r, p_val = compute_correlation(csa, edss)
            print(f'{site}, {phenotype}: Correlation EDSS vs CSA: r={r}, p-value{format_pvalue(p_val, alpha=0.05)}')
            # Plot individual scatter plots
            ax[index].scatter(csa, edss, color=color_pallete[color], alpha=.8, label=phenotype, s=100)
            x_vals, y_vals = compute_regression(csa, edss)
            ax[index].plot(x_vals, y_vals, '--', color=color_pallete[color], alpha=.8, linewidth=3)
            if site == 'all':
                ax[index].set_title(site_to_vendor_title[site], fontsize=FONTSIZE_CORR, fontweight='bold')
            else:
                ax[index].set_title(site_to_vendor_title[site], fontsize=FONTSIZE_CORR)
            if index > 2:
                ax[index].set_xlabel('CSA [$mm^2$]', fontsize=FONTSIZE_CORR)

            # # Set fixed number of y-ticks
            # xmin, xmax = ax[index].get_xlim()
            # custom_ticks = np.linspace(50, xmax, 5, dtype=int)
            # ax[index].set_xticks(custom_ticks)
            # ax[index].set_xticklabels(custom_ticks)

            if index == 0 or index == 3:
                ax[index].set_ylabel('EDSS', fontsize=FONTSIZE_CORR)
            # Increase size of xticks and yticks
            plt.setp(ax[index].xaxis.get_majorticklabels(), fontsize=FONTSIZE_CORR)
            plt.setp(ax[index].yaxis.get_majorticklabels(), fontsize=FONTSIZE_CORR)
            # Show legend only for the last axis
            # loc is used to move the legend below the text with corr coef and pval
            if index == 5:
                #legend = ax[index].legend(fontsize=FONTSIZE_CORR - 5, loc=(.65, .63))      # 5 sites
                legend = ax[index].legend(fontsize=FONTSIZE_CORR - 5, loc=(.68, .66))      # all + 5 sites
                #legend = ax[index].legend(fontsize=FONTSIZE_CORR - 5, loc=(.68, .75))       # all + 5 sites; no text
                frame = legend.get_frame()  # sets up for color, edge, and transparency
                frame.set_edgecolor('black')  # edge color of legend
    plt.subplots_adjust(wspace=-0.1)
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
        if site == 'all':
            continue        # TODO: include 'all' to the df
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
    canproco_pd = pd.merge(canproco_pd, canproco_participants_pd[['participant_id', 'pathology', 'phenotype', 'edss_M0']],
                           how='left', left_on='subject_id', right_on='participant_id')

    # get MS patients with 'n/a' for EDSS
    # canproco_pd[(canproco_pd['edss_M0'].isna()) & (canproco_pd['pathology'] == 'MS')]['subject_id']

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

    # Compute and plot correlation between EDSS and CSA persite
    fname_fig = args.i_canproco.replace('.csv', '_correlation_persite.png')
    create_csa_edss_correlation_figure_persite(canproco_pd, fname_fig)


if __name__ == "__main__":
    main()
