"""
In this file we analyse the results of the segmentation of the MS lesion on the spinal cord.
The objective is to output the number of lesions segmented on the spinal cord for each patient.
Also, we want to output the segmentation file of the lesions in different color

Usage:
    python lesion_seg_analysis.py -i <input_image> -seg <segmentation> -o <output_folder>

Args:
    -i/--input_image: path to the image
    -seg/--segmentation: path to the segmentation
    -o/--output_folder: path to the output folder
    --plot: whether to plot the results or not

Returns:
    None

Todo:
    *

Pierre-Louis Benveniste
"""

import os
import argparse
from pathlib import Path
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN


def get_parser():
    """
    This function parses the arguments given to the script.

    Args:
        None

    Returns:
        parser: parser containing the arguments
    """

    parser = argparse.ArgumentParser(description='Analyse the results of the segmentation of the MS lesion on the spinal cord.')
    parser.add_argument('-i', '--input_image', required=True,
                        help='Path to the image')
    parser.add_argument('-seg', '--segmentation', required=True,
                        help='Path to the segmentation')
    parser.add_argument('-o', '--output_folder', required=True,
                        help='Path to the output folder')
    parser.add_argument('--plot', action='store_true',
                        help='Whether to plot the results or not')
    return parser


def main():
    """
    This function is the main function of the script.

    Args:
        None
    
    Returns:
        None
    """
    #get the parser
    parser = get_parser()
    args = parser.parse_args()

    #load image
    img = nib.load(args.input_image)
    img_data = img.get_fdata()

    #load segmentation
    seg = nib.load(args.segmentation)
    seg_data = seg.get_fdata()

    #perform clustering on the entire segmentation volume
    ##first we modify the seg data
    X = []
    Y = []
    Z = []
    for x in range(seg_data.shape[0]):
        for y in range(seg_data.shape[1]):
            for z in range(seg_data.shape[2]):
                if seg_data[x,y,z] != 0:
                    X.append(x)
                    Y.append(y)
                    Z.append(z)
    coords = np.stack((X,Y,Z), axis=1)

    ##then we perform the clustering using DBSCAN
    db = DBSCAN(eps=10, min_samples=5).fit(coords)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

    print('Estimated number of clusters: %d' % n_clusters_)

    plot = args.plot
    if plot:
        #build color dictionnary
        colors = {}
        for i in range(n_clusters_):
            colors[i] = np.random.rand(3,)
        colors[-1] = [0,0,0]

        #plot the clusters for each slice
        fig, axs = plt.subplots(ncols=seg_data.shape[2], nrows=1,figsize=(20,3))
        for i in range(seg_data.shape[2]):
            slice_coords = coords[coords[:,2] == i]
            slice_labels = labels[coords[:,2] == i]
            axs[i].scatter(slice_coords[:,0], slice_coords[:,1], color=[colors[x] for x in slice_labels])
            # plt.savefig(os.path.join(args.output_folder, f'cluster_slice_{i}.png'))
            # plt.close()
            axs[i].set_xlim(0,seg_data.shape[0])   
            axs[i].set_ylim(0,seg_data.shape[1])
        plt.show()

    #saving the results in a nifti file
    #first we create the new segmentation
    new_seg_data = np.zeros_like(seg_data)
    for i in range(len(labels)):
        new_seg_data[coords[i,0], coords[i,1], coords[i,2]] = labels[i] + 1
    #then we save it
    new_seg = nib.Nifti1Image(new_seg_data, seg.affine, seg.header)
    nib.save(new_seg, os.path.join(args.output_folder, 'clustered_seg.nii.gz'))

    #for each lesion calculate volume and get center
    ##first we get the volume of one voxel
    voxel_volume = np.prod(seg.header.get_zooms())
    ##then we get the volume of each lesion and its center
    lesion_volumes = []
    lesion_centers = []
    for i in range(n_clusters_):
        lesion_volumes.append(len(labels[labels == i])*voxel_volume)
        lesion_centers.append(np.mean(coords[labels == i], axis=0))
    
    #save the results in a text file
    with open(os.path.join(args.output_folder, 'lesion_analysis.txt'), 'w') as f:
        f.write(f'Number of lesions: {n_clusters_}\n')
        f.write('Volume and center of each lesion (mm3):\n')
        for i in range(n_clusters_):
            f.write(f'Lesion {i+1} : volume: {round(lesion_volumes[i],2)} mm3, center: {lesion_centers[i]}\n')
    return None


if __name__ == '__main__':
    main()

