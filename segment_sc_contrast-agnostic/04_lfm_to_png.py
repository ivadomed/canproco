"""
Generate axial png images from the Lesion Frequency Maps (LFM) for each vertebral level.

The script requires the SCT conda environment to be activated:
    source ${SCT_DIR}/python/etc/profile.d/conda.sh
    conda activate venv_sct

Author: Jan Valosek
Inspired by https://github.com/neuropoly/lesion-mapping/blob/master/lfm_2_png.py from Charley Gros
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import find_contours
from scipy.ndimage.morphology import binary_dilation

from spinalcordtoolbox.image import Image


def get_parser():
    """
    parser function
    """

    parser = argparse.ArgumentParser(
        description='Generate axial png images from the Lesion Frequency Maps (LFM).',
        prog=os.path.basename(__file__).strip('.py')
    )
    parser.add_argument(
        '-lfm-path',
        metavar="<file>",
        required=True,
        type=str,
        help='Path to the LFM file. Example: /results/spinalcord_LFM_MS.nii.gz'
    )
    parser.add_argument(
        '-thr',
        metavar="<file>",
        required=False,
        type=str,
        help='Threshold for the LFM in range of 0-1. Default: 0.15 (meaning 15%)',
        default=0.1
    )
    parser.add_argument(
        '-ofolder',
        metavar="<folder>",
        required=True,
        type=str,
        help='Path to the output folder. Example: /results/'
    )

    return parser


def load_data(path, thr_bin):
    """
    Load data from the path and threshold it.
    """
    img = Image(path).change_orientation('RPI')
    # Threshold
    if thr_bin > 0:
        data = img.data
        data = (data > thr_bin).astype(np.int)
    else:
        data = (img.data > 0).astype(np.int)
    del img
    return data


def rescale_rot(img, rescale):
    """
    Rescale and rotate the image.
    """
    img = np.repeat(img, rescale, axis=0)
    img = np.repeat(img, rescale, axis=1)
    img = np.rot90(img)
    return img


def combine_img_w_bkg(img, bkg, gm, rescale, thr, fname_out, linewidth=4):
    """
    Combine the LFM with the background image (PAM50_t2) and the grey matter segmentation.
    """
    i_zero, i_nonzero = np.where(img == 0.0), np.nonzero(img)

    img_jet = plt.cm.jet(plt.Normalize(vmin=0, vmax=thr)(img))
    img_jet[i_zero] = 0.0

    bkg_grey = plt.cm.binary_r(plt.Normalize(vmin=np.amin(bkg), vmax=np.amax(bkg))(bkg))

    img_out = np.copy(bkg_grey)
    img_out[i_nonzero] = img_jet[i_nonzero]

    img_out = rescale_rot(img_out, rescale)
    gm = rescale_rot(gm, rescale)

    ratio_shape = img_out.shape[0] * 1. / img_out.shape[1]
    plt.figure(figsize=(10, 10 * ratio_shape))
    plt.subplot(1, 1, 1)
    plt.axis("off")
    plt.imshow(img_out, interpolation='nearest', aspect='auto')

    #print(np.unique(gm))
    gm_dilated = binary_dilation(gm)
    #print(np.unique(gm))
    contours = find_contours(gm_dilated, .5)
    for n, contour in enumerate(contours):
        plt.plot(contour[:, 1], contour[:, 0], 'white', linewidth=linewidth)
    for n, contour in enumerate(contours):
        plt.plot(contour[:, 1], contour[:, 0], 'black', linewidth=linewidth // 2)

    plt.savefig(fname_out, dpi=800)
    plt.close()
    print('Saved: ' + fname_out)


def save_colormap(fname_out, cmap='jet'):
    """
    Save matplotlib colormap as png
    """
    fig = plt.figure(figsize=[10, 1])
    ax = fig.add_subplot(111)

    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))

    ax.imshow(gradient, aspect='auto', cmap=cmap)

    ax.set_axis_off()
    fig.savefig(fname_out)
    plt.close()


def load_PAM50_gm():
    """
    Load PAM50 grey matter mask
    """
    path_template = os.path.join(os.environ.get('SCT_DIR'), 'data', 'PAM50', 'template')
    data = load_data(os.path.join(path_template, 'PAM50_gm.nii.gz'), thr_bin=0.7)

    return data


def main():
    # Parse the command line arguments
    parser = get_parser()
    args = parser.parse_args()

    lfm_path = args.lfm_path
    thr = float(args.thr)
    ofolder = args.ofolder

    if not os.path.isdir(ofolder):
        os.makedirs(ofolder)

    img_img = Image(lfm_path).change_orientation('RPI')
    img = img_img.data
    del img_img

    path_pam50_t2 = os.path.join(os.environ.get('SCT_DIR'), 'data', 'PAM50', 'template', 'PAM50_t2.nii.gz')
    img_pam50_t2 = Image(path_pam50_t2)
    backgroud = img_pam50_t2.data
    del img_pam50_t2
    x_shape, y_shape, z_shape = backgroud.shape
    x_mean, y_mean = x_shape // 2, y_shape // 2
    # PAM50_t2
    backgroud = backgroud[x_mean - 25:x_mean + 25, y_mean - 25:y_mean + 25, :]
    # LFM
    img = img[x_mean - 25:x_mean + 25, y_mean - 25:y_mean + 25, :]

    gm_mask = load_PAM50_gm()
    gm_mask = gm_mask[x_mean - 25:x_mean + 25, y_mean - 25:y_mean + 25, :]

    path_pam50_disc = os.path.join(os.environ.get('SCT_DIR'), 'data', 'PAM50', 'template', 'PAM50_label_disc.nii.gz')
    img_lvl = Image(path_pam50_disc)
    data_lvl = img_lvl.data
    del img_lvl

    lvl_z_lst = [np.where(data_lvl == lvl)[2][0] for lvl in np.unique(data_lvl) if lvl in range(1, 9)]

    img_lst, bkg_lst, gm_lst = [], [], []
    # Loop across vertebral levels
    for lvl_idx in range(len(lvl_z_lst) - 1):
        # Get min top and bottom slices for current level
        z_bot, z_top = lvl_z_lst[lvl_idx + 1], lvl_z_lst[lvl_idx] + 1
        # Compute middle slice for current level
        z_mid_lvl = (z_top + z_bot) // 2
        # Average LFM across slices for current level
        img_lst.append(np.mean(img[:, :, z_bot:z_top], axis=2))
        # Get middle slice for background (PAM50_t2)
        bkg_lst.append(backgroud[:, :, z_mid_lvl])
        gm_lst.append(gm_mask[:, :, z_mid_lvl])
    # Construct list of of vert levels C1-C7
    pref_lst = ['C' + str(i) for i in range(1, 8)]

    for img_cur, bkg_cur, gm_cur, pref_cur in zip(img_lst, bkg_lst, gm_lst, pref_lst):
        fname_out_cur = os.path.join(ofolder, lfm_path.split('/')[-1].replace('.nii.gz', '_' + pref_cur + '.png'))
        combine_img_w_bkg(img_cur, bkg_cur, gm_cur, rescale=4, thr=thr, fname_out=fname_out_cur)

    # Save colormap
    save_colormap(os.path.join(ofolder, 'jet_0_' + str(int(thr * 100)) + '.png'))


if __name__ == "__main__":
    main()
