"""
This python file is used to analyse the evolution of the lesions over time.
It first performs a registration of the images and their segmentation on the two time points.
Then, it performs clustering of the lesions on each image and computes the volume and the center of each lesion. 
Finally, it compares the lesions between the two time points and computes the evolution of the lesions.

Args:
    -i-1, --input-first_image: path to the first image 
    -i-2, --input-second_image: path to the second image
    -seg-1, --segmentation-first_image: path to the lesion segmentation of the first image
    -seg-2, --segmentation-second_image: path to the lesion segmentation of the second image
    -o, --output_folder: path to the output folder
    --plot: whether to plot the results or not

Returns:
    None

Example:
    python time_point_lesion_evolution.py -i-1 path/to/first/image -i-2 path/to/second/image -seg-1 path/to/first/segmentation -seg-2 path/to/second/segmentation -o path/to/output/folder
    python3 lesion_analysis/time_point_lesion_evolution.py -i-1 /Users/plbenveniste/Desktop/lesion_comparison_copy/sub-cal072_ses-M0_STIR.nii.gz -i-2 /Users/plbenveniste/Desktop/lesion_comparison_copy/sub-cal072_ses-M12_STIR.nii.gz -seg-1 /Users/plbenveniste/Desktop/lesion_comparison_copy/sub-cal072_ses-M0_STIR_lesion-manual.nii.gz -seg-2 /Users/plbenveniste/Desktop/lesion_comparison_copy/M12_inference_results.nii.gz -o /Users/plbenveniste/Desktop/lesion_comparison_copy/output

To do:
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
import networkx as nx
from scipy.optimize import linear_sum_assignment

def get_parser():
    """
    This function parses the arguments given to the script.

    Args:
        None

    Returns:
        parser: parser containing the arguments
    """

    parser = argparse.ArgumentParser(description='Analyse the results of the segmentation of the MS lesion on the spinal cord.')
    parser.add_argument('-i-1', '--input-first_image', required=True,
                        help='Path to the first image')
    parser.add_argument('-i-2', '--input-second_image', required=True,
                        help='Path to the second image')
    parser.add_argument('-seg-1', '--segmentation-first_image', required=True,
                        help='Path to the lesion segmentation of the first image')
    parser.add_argument('-seg-2', '--segmentation-second_image', required=True,
                        help='Path to the lesion segmentation of the second image')
    parser.add_argument('-o', '--output_folder', required=True,
                        help='Path to the output folder')
    parser.add_argument('--show', action='store_true',
                        help='Whether to show the results or not')
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

    #get the images and the segmentations
    first_image = nib.load(args.input_first_image)
    second_image = nib.load(args.input_second_image)
    first_lesion_segmentation = nib.load(args.segmentation_first_image)
    second_lesion_segmentation = nib.load(args.segmentation_second_image)

    """
    #first we segment the spinal cord on the two images using sct_deepseg_sc
    os.system('sct_deepseg_sc -i ' + args.input_first_image + ' -c t2 -o ' + args.output_folder + '/first_image_sc_segmentation.nii.gz')
    os.system('sct_deepseg_sc -i ' + args.input_second_image + ' -c t2 -o' + args.output_folder + '/second_image_sc_segmentation.nii.gz')

    #then we compute the vertebral levels of the images using sct_label_vertebrae
    os.system('sct_label_vertebrae -i ' + args.input_first_image + ' -s ' + args.output_folder + '/first_image_sc_segmentation.nii.gz' + ' -c t2 -ofolder ' + args.output_folder)
    os.system('sct_label_vertebrae -i ' + args.input_second_image + ' -s ' + args.output_folder + '/second_image_sc_segmentation.nii.gz' + ' -c t2 -ofolder ' + args.output_folder)

    first_label = nib.load(args.output_folder + '/first_image_sc_segmentation_labeled.nii.gz')
    second_label = nib.load(args.output_folder + '/second_image_sc_segmentation_labeled.nii.gz')
    first_label_data = first_label.get_fdata()
    second_label_data = second_label.get_fdata()
    common_labels = np.intersect1d(np.unique(first_label_data), np.unique(second_label_data))

    #we remove labels that are not in common in both images
    for label in np.unique(first_label_data):
        if label not in common_labels:
            first_label_data[first_label_data == label] = 0
    #we save the new segmentation
    first_label = nib.Nifti1Image(first_label_data, first_label.affine, first_label.header)
    nib.save(first_label, args.output_folder + '/first_image_sc_segmentation_labeled.nii.gz')
    #same for the second image
    for label in np.unique(second_label_data):
        if label not in common_labels:
            second_label_data[second_label_data == label] = 0
    second_label = nib.Nifti1Image(second_label_data, second_label.affine, second_label.header)
    nib.save(second_label, args.output_folder + '/second_image_sc_segmentation_labeled.nii.gz')
    
    #then we register the images and the segmentations with the first image as reference using sct_register_multimodal using both the spinal cord segmentation and the vertebral levels
    parameters = 'step=0,type=label,dof=Tx_Ty_Tz:step=1,type=im,algo=dl'
    os.system('sct_register_multimodal -i ' + args.input_second_image + ' -d ' + args.input_first_image + ' -ilabel ' + args.output_folder + '/second_image_sc_segmentation_labeled.nii.gz' 
                + ' -dlabel ' + args.output_folder + '/first_image_sc_segmentation_labeled.nii.gz' + ' -ofolder ' + args.output_folder + ' -owarp ' + args.output_folder + '/warp_M12_to_M0.nii.gz'
                 + ' -owarpinv ' + args.output_folder + '/warp_M0_to_M12.nii.gz' + ' -param ' + parameters + ' -x linear ')

    #we apply the warping to the segmentation of the second image
    os.system('sct_apply_transfo -i ' + args.segmentation_second_image + ' -d ' + args.segmentation_first_image + ' -w ' + args.output_folder + '/warp_M12_to_M0.nii.gz' + ' -o ' + args.output_folder + '/second_lesion_seg_registered.nii.gz' + ' -x linear')
    """

    #we first look at the values in the reg lesion  segmentation
    second_lesion_seg = nib.load(args.output_folder + '/second_lesion_seg_registered.nii.gz')
    second_lesion_seg_data = second_lesion_seg.get_fdata()
    #we replace all values above a threshold by 1
    second_lesion_seg_data[second_lesion_seg_data > 1e-5] = 1

    #we then load data from the first lesion segmentation
    first_lesion_seg = nib.load(args.segmentation_first_image)
    first_lesion_seg_data = first_lesion_seg.get_fdata()

    # we then perform clustering on the two images
    # we first get the coordinates of the lesions
    first_lesion_seg_coordinates = np.argwhere(first_lesion_seg_data == 1)
    second_lesion_seg_coordinates = np.argwhere(second_lesion_seg_data == 1)
    # we then perform clustering on the coordinates
    clustering = DBSCAN(eps=10, min_samples=5).fit(first_lesion_seg_coordinates)
    first_lesion_seg_labels = clustering.labels_
    clustering = DBSCAN(eps=10, min_samples=5).fit(second_lesion_seg_coordinates)
    second_lesion_seg_labels = clustering.labels_

    #output the nifti images of the lesions where each lesion has a different value
    first_lesion_seg_data2 = np.zeros(first_lesion_seg_data.shape)
    second_lesion_seg_data2 = np.zeros(second_lesion_seg_data.shape)
    for i,voxel in enumerate(first_lesion_seg_coordinates):
        first_lesion_seg_data2[voxel[0], voxel[1], voxel[2]] = first_lesion_seg_labels[i]+1
    for i,voxel in enumerate(second_lesion_seg_coordinates):
        second_lesion_seg_data2[voxel[0], voxel[1], voxel[2]] = second_lesion_seg_labels[i]+1
    first_lesion_seg2 = nib.Nifti1Image(first_lesion_seg_data2, first_lesion_seg.affine, first_lesion_seg.header)
    second_lesion_seg2 = nib.Nifti1Image(second_lesion_seg_data2, second_lesion_seg.affine, second_lesion_seg.header)
    nib.save(first_lesion_seg2, args.output_folder + '/first_lesion_seg_clustered.nii.gz')
    nib.save(second_lesion_seg2, args.output_folder + '/second_lesion_seg_clustered.nii.gz')

    #for each lesion of the second time point we consider that it corresponds to the lesion of the first time point with which it overlaps more than 50%
    matching_table = np.empty(len(np.unique(second_lesion_seg_data2))-1)
    for lesion_2nd in np.unique(second_lesion_seg_data2):
        if lesion_2nd != 0:
            for lesion_1st in np.unique(first_lesion_seg_data2):
                if lesion_1st !=0:
                    #we get the coordinates of the lesion
                    second_lesion_seg_coordinates2 = np.argwhere(second_lesion_seg_data2 == lesion_2nd)
                    #we get the coordinates of the lesion in the first time point
                    first_lesion_seg_coordinates = np.argwhere(first_lesion_seg_data2 == lesion_1st)
                    #we compute the intersection of the two sets of coordinates
                    second_set = set( tuple(points) for points in second_lesion_seg_coordinates2)
                    first_set = set(tuple(points) for points in  first_lesion_seg_coordinates)
                    intersection = len(second_set.intersection(first_set))
                    ratio = intersection/len(first_set)
                    if ratio > 0.5:
                        matching_table[int(lesion_2nd-1)] = int(lesion_1st)
                    else:
                        matching_table[int(lesion_2nd-1)] = None

    #we then compute the center of each lesion
    center_2nd = []
    for lesion_2nd in np.unique(second_lesion_seg_data2):
        if lesion_2nd != 0:
            second_lesion_seg_coordinates2 = np.argwhere(second_lesion_seg_data2 == lesion_2nd)
            center_2nd.append(np.mean(second_lesion_seg_coordinates2, axis=0))
    center_2nd = np.array(center_2nd)
    center_1st = []
    for lesion_1st in np.unique(first_lesion_seg_data2):
        if lesion_1st != 0:
            first_lesion_seg_coordinates = np.argwhere(first_lesion_seg_data2 == lesion_1st)
            center_1st.append(np.mean(first_lesion_seg_coordinates, axis=0))
    center_1st = np.array(center_1st)
 
    #we then perfrom lesion clustering on the original seg file from the second time point
    real_second_lesion_seg = nib.load(args.segmentation_second_image)
    real_second_lesion_seg_data = real_second_lesion_seg.get_fdata()
    real_second_lesion_seg_coordinates = np.argwhere(real_second_lesion_seg_data == 1)
    clustering = DBSCAN(eps=10, min_samples=5).fit(real_second_lesion_seg_coordinates)
    real_second_lesion_seg_labels = clustering.labels_
    #for these we conpute the center
    real_second_lesion_seg_data2 = np.zeros(real_second_lesion_seg_data.shape)
    for i,voxel in enumerate(real_second_lesion_seg_coordinates):
        real_second_lesion_seg_data2[voxel[0], voxel[1], voxel[2]] = real_second_lesion_seg_labels[i]+1
    center_real_2nd = []
    for lesion_2nd in np.unique(real_second_lesion_seg_data2):
        if lesion_2nd != 0:
            real_second_lesion_seg_coordinates2 = np.argwhere(real_second_lesion_seg_data2 == lesion_2nd)
            center_real_2nd.append(np.mean(real_second_lesion_seg_coordinates2, axis=0))
    center_real_2nd = np.array(center_real_2nd)
    
    #we then perform matching between lesions of the second time point warpped and real lesions of the secont time point
    matching_table2 = np.empty(len(np.unique(real_second_lesion_seg_data2))-1)

    # Function to compute distance between two lesions (you can customize this based on your features)
    def lesion_distance(lesion1, lesion2):
        return np.sqrt((lesion1[0] - lesion2[0])**2 + (lesion1[1] - lesion2[1])**2 + (lesion1[2] - lesion2[2])**2)

    # Create a graph for each time point
    graph_warp = nx.Graph()
    graph_t2 = nx.Graph()

    # Add nodes (lesions) to each graph
    graph_warp.add_nodes_from(range(len(center_2nd)))
    graph_t2.add_nodes_from(range(len(center_real_2nd)))

    # Add edges between nodes with similarity based on lesion distance
    for i, lesion1 in enumerate(center_2nd):
        for j, lesion2 in enumerate(center_real_2nd):
            distance = lesion_distance(lesion1, lesion2)
            graph_warp.add_edge(i, j, weight=distance)
            graph_t2.add_edge(j, i, weight=distance)  # Adding bidirectional edges

    # Use Hungarian algorithm for graph matching
    cost_matrix = np.array(nx.to_numpy_array(graph_warp))
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # Create temporal links based on matching results
    warping_links = [(i, j) for i, j in zip(range(len(center_2nd)), col_ind) if cost_matrix[i, j] < np.inf]

    print("Temporal Links:",warping_links)

if __name__ == '__main__':
    main()