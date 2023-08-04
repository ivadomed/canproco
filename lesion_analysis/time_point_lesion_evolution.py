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


def perform_registration(path_first_image, path_second_image, path_first_lesion_seg, path_second_lesion_seg, path_output_folder):
    """
    This function performs the registration of the images and the segmentations using spinal cord toolbox

    Args:
        path_first_image: path to the first image
        path_second_image: path to the second image
        path_first_lesion_seg: path to the lesion segmentation of the first image
        path_second_lesion_seg: path to the lesion segmentation of the second image
        path_output_folder: path to the output folder
    
    Returns:
        path_warp: path to the warp file
        path_warpinv: path to the inverse warp file
        path_second_image_registered: path to the second image registered
        path_second_lesion_seg_registered: path to the second lesion segmentation registered
    """

    #first we segment the spinal cord on the two images using sct_deepseg_sc
    os.system('sct_deepseg_sc -i ' + path_first_image + ' -c t2 -o ' + path_output_folder + '/first_image_sc_seg.nii.gz')
    os.system('sct_deepseg_sc -i ' + path_second_image + ' -c t2 -o' + path_output_folder + '/second_image_sc_seg.nii.gz')

    #then we compute the vertebral levels of the images using sct_label_vertebrae
    os.system('sct_label_vertebrae -i ' + path_first_image + ' -s ' + path_output_folder + '/first_image_sc_seg.nii.gz' + ' -c t2 -ofolder ' + path_output_folder)
    os.system('sct_label_vertebrae -i ' + path_second_image + ' -s ' + path_output_folder + '/second_image_sc_seg.nii.gz' + ' -c t2 -ofolder ' + path_output_folder)

    #we remove labels that are not in common in both images
    ##we first get the labels of the two images
    first_label = nib.load(path_output_folder + '/first_image_sc_seg_labeled.nii.gz')
    second_label = nib.load(path_output_folder + '/second_image_sc_seg_labeled.nii.gz')
    first_label_data = first_label.get_fdata()
    second_label_data = second_label.get_fdata()
    common_labels = np.intersect1d(np.unique(first_label_data), np.unique(second_label_data))

    #we then remove the labels that are not in common
    for label in np.unique(first_label_data):
        if label not in common_labels:
            first_label_data[first_label_data == label] = 0
    #we overwrite the segmentation to keep only common labels
    first_label = nib.Nifti1Image(first_label_data, first_label.affine, first_label.header)
    nib.save(first_label, path_output_folder + '/first_image_sc_seg_labeled.nii.gz')
    #same for the second image
    for label in np.unique(second_label_data):
        if label not in common_labels:
            second_label_data[second_label_data == label] = 0
    second_label = nib.Nifti1Image(second_label_data, second_label.affine, second_label.header)
    nib.save(second_label, path_output_folder + '/second_image_sc_seg_labeled.nii.gz')
    
    #then we register the images and the segmentations with the first image as reference using sct_register_multimodal using both the spinal cord segmentation and the vertebral levels
    parameters = 'step=0,type=label,dof=Tx_Ty_Tz:step=1,type=im,algo=dl'
    os.system('sct_register_multimodal -i ' + path_second_image + ' -d ' + path_first_image + ' -ilabel ' + path_output_folder + '/second_image_sc_seg_labeled.nii.gz' 
                + ' -dlabel ' + path_output_folder + '/first_image_sc_seg_labeled.nii.gz' + ' -ofolder ' + path_output_folder + ' -owarp ' + path_output_folder + '/warp_2nd_to_1st.nii.gz'
                 + ' -owarpinv ' + path_output_folder + '/warp_1st_to_2nd.nii.gz' + ' -param ' + parameters + ' -x linear ')

    #we apply the warping to the segmentation of the second image
    os.system('sct_apply_transfo -i ' + path_second_lesion_seg + ' -d ' + path_first_lesion_seg + ' -w ' + path_output_folder + '/warp_2nd_to_1st.nii.gz' 
              + ' -o ' + path_output_folder + '/second_lesion_seg_registered.nii.gz' + ' -x linear')

    return path_output_folder + '/warp_2nd_to_1st.nii.gz', path_output_folder + '/warp_1st_to_2nd.nii.gz', path_output_folder + '/second_image_sc_seg_registered.nii.gz',path_output_folder + '/second_lesion_seg_registered.nii.gz'


def temporal_lesion_matching_on_reg_image(path_first_lesion_seg, path_second_lesion_seg_reg):
    """
    This function performs temporal lesion matching using the 1st time point lesion segmentation and the 2nd time point lesion segmentation registered to the 1st time point.

    Args:
        path_first_lesion_seg: path to the lesion segmentation of the first image
        path_second_lesion_seg_reg: path to the lesion segmentation of the second image registered to the first image

    Returns:
        matching_table: table containing the matching between lesions of the first time point and lesions of the second time point
        center_1st: coordinates of the center of the lesions of the first time point
        center_2nd: coordinates of the center of the lesions of the second time point
        1st_data_labeled: segmentation of the first time point with each lesion in a different color
        2nd_data_labeled: segmentation of the second time point with each lesion in a different color
    """

    #because registration modifies the values of the segmentation, we first replace all values above a threshold by 1
    second_lesion_seg_reg = nib.load(path_second_lesion_seg_reg)
    second_lesion_seg_reg_data = second_lesion_seg_reg.get_fdata()
    second_lesion_seg_reg_data[second_lesion_seg_reg_data > 1e-5] = 1

    #we then load data from the first lesion segmentation
    first_lesion_seg = nib.load(path_first_lesion_seg)
    first_lesion_seg_data = first_lesion_seg.get_fdata()

    # we then perform clustering on the two images
    # we first get the coordinates of the lesions
    first_lesion_seg_coordinates = np.argwhere(first_lesion_seg_data == 1)
    second_lesion_seg_reg_coordinates = np.argwhere(second_lesion_seg_reg_data == 1)
    # we then perform clustering on the coordinates
    clustering_1st = DBSCAN(eps=10, min_samples=5).fit(first_lesion_seg_coordinates)
    first_lesion_seg_labels = clustering_1st.labels_
    clustering_2nd = DBSCAN(eps=10, min_samples=5).fit(second_lesion_seg_reg_coordinates)
    second_lesion_seg_reg_labels = clustering_2nd.labels_

    #create the same tables with the lesions in different colors (voxel value = lesion ID)
    first_lesion_seg_data2 = np.zeros(first_lesion_seg_data.shape)
    second_lesion_seg_reg_data2 = np.zeros(second_lesion_seg_reg_data.shape)
    for i,voxel in enumerate(first_lesion_seg_coordinates):
        first_lesion_seg_data2[voxel[0], voxel[1], voxel[2]] = first_lesion_seg_labels[i]+1
    for i,voxel in enumerate(second_lesion_seg_reg_coordinates):
        second_lesion_seg_reg_data2[voxel[0], voxel[1], voxel[2]] = second_lesion_seg_reg_labels[i]+1

    #for each lesion of the second time point we consider that it corresponds to the lesion of the first time point with which it overlaps more than 50%
    matching_table = {}
    for lesion_2nd in np.unique(second_lesion_seg_reg_data2):
        if lesion_2nd != 0:
            for lesion_1st in np.unique(first_lesion_seg_data2):
                if lesion_1st !=0:
                    #we get the coordinates of the lesion
                    second_lesion_seg_coordinates2 = np.argwhere(second_lesion_seg_reg_data2 == lesion_2nd)
                    #we get the coordinates of the lesion in the first time point
                    first_lesion_seg_coordinates = np.argwhere(first_lesion_seg_data2 == lesion_1st)
                    #we compute the overlap of the two sets of coordinates
                    second_set = set( tuple(points) for points in second_lesion_seg_coordinates2)
                    first_set = set(tuple(points) for points in  first_lesion_seg_coordinates)
                    intersection = len(second_set.intersection(first_set))
                    ratio = intersection/len(first_set)
                    matching_table[int(lesion_2nd)] = {}
                    if ratio > 0.5:
                        matching_table[int(lesion_2nd)]["Lesion_1st_timepoint"] = int(lesion_1st)
                        matching_table[int(lesion_2nd)]["Overlapping_ratio"] = round(ratio,2)
                    else:
                        matching_table[int(lesion_2nd)]["Lesion_1st_timepoint"] = None
                        matching_table[int(lesion_2nd)]["Overlapping_ratio"] = None

    #we then compute the center of each lesion
    center_2nd = {}
    for lesion_2nd in np.unique(second_lesion_seg_reg_data2):
        if lesion_2nd != 0:
            second_lesion_seg_reg_coordinates2 = np.argwhere(second_lesion_seg_reg_data2 == lesion_2nd)
            center_2nd[int(lesion_2nd)] = np.mean(second_lesion_seg_reg_coordinates2, axis=0)
    center_1st = {}
    for lesion_1st in np.unique(first_lesion_seg_data2):
        if lesion_1st != 0:
            first_lesion_seg_coordinates = np.argwhere(first_lesion_seg_data2 == lesion_1st)
            center_1st[int(lesion_1st)] = np.mean(first_lesion_seg_coordinates, axis=0)

    #we first overwrite the segmentation files with the lesions in different colors
    final_first_lesion_seg = np.zeros(first_lesion_seg_data2.shape)
    final_second_lesion_seg_reg = np.zeros(second_lesion_seg_reg_data2.shape)
    for i,voxel in enumerate(first_lesion_seg_coordinates):
        final_first_lesion_seg[voxel[0], voxel[1], voxel[2]] = matching_table[first_lesion_seg_data2[voxel[0], voxel[1], voxel[2]]]["Lesion_1st_timepoint"]

    for i,lesion in enumerate(matching_table):
        print(i)
        print(lesion)
        print(matching_table[lesion])
        first_lesion_seg_data2.replace(matching_table[lesion]["Lesion_1st_timepoint"], lesion)
    #print(first_lesion_seg_data2)


    #first_lesion_seg2 = nib.Nifti1Image(first_lesion_seg_data2, first_lesion_seg.affine, first_lesion_seg.header)
    #second_lesion_seg2 = nib.Nifti1Image(second_lesion_seg_data2, second_lesion_seg.affine, second_lesion_seg.header)
    #nib.save(first_lesion_seg2, args.output_folder + '/first_lesion_seg_clustered.nii.gz')
    #nib.save(second_lesion_seg2, args.output_folder + '/second_lesion_seg_clustered.nii.gz')

    return matching_table, center_1st, center_2nd, first_lesion_seg_data2, second_lesion_seg_reg_data2


def lesion_distance(lesion1, lesion2, volume1, volume2):
    """
    This function computes the distance betweeen two lesions for the graph matching.

    Args:
        lesion1: coordinates of the first lesion
        lesion2: coordinates of the second lesion
    
    Returns:
        distance: distance between the two lesions
    """
    pos_distance = np.sqrt((lesion1[0] - lesion2[0])**2 + (lesion1[1] - lesion2[1])**2 + (lesion1[2] - lesion2[2])**2)
    vol_distance = np.abs(volume1 - volume2)

    return  pos_distance + vol_distance


def volume_lesions(labeled_data, voxel_volume):
    """
    This function computes the volume of the lesions.

    Args:
        labeled_data: data with each lesion in a different color
        voxel_volume: volume of a voxel
    
    Returns:
        volumes: volume of the lesions
    """

    volumes = []
    for lesion in np.unique(labeled_data):
        if lesion != 0:
            volume = len(np.argwhere(labeled_data == lesion))*voxel_volume
            volumes.append((lesion,volume))
    return np.array(volumes)  


def graph_matching(center_2nd_orig, center_2nd_reg, second_orig_volumes, second_reg_volumes):
    """
    This function performs graph matching between the lesions of the second time point registered and the lesions of the second time point.
    The node distances take into account the distance between the lesions and the volume of the lesions.

    Args:
        center_2nd: coordinates of the center of the lesions of the second time point
        center_2nd_reg: coordinates of the center of the lesions of the second time point registered
        second_orig_volumes: volume of the lesions of the second time point
        second_reg_volumes: volume of the lesions of the second time point registered
    
    Returns:
        matching_table: table containing the matching between lesions of the second time point and lesions of the second time point registered
    """

    # Create a graph for each time point
    graph_t2_orig = nx.Graph()
    graph_t2_reg = nx.Graph()

    # Add nodes (lesions) to each graph
    graph_t2_orig.add_nodes_from(range(len(center_2nd_orig)))
    graph_t2_reg.add_nodes_from(range(len(center_2nd_reg)))

    # Add edges between nodes with similarity based on lesion distance
    for i, lesion_reg in enumerate(center_2nd_reg):
        for j, lesion_orig in enumerate(center_2nd_orig):
            distance = lesion_distance(center_2nd_reg[lesion_reg], center_2nd_orig[lesion_orig], second_reg_volumes[i][1], second_orig_volumes[j][1])
            graph_t2_reg.add_edge(i, j, weight=distance)
            graph_t2_orig.add_edge(j, i, weight=distance)  # Adding bidirectional edges

    # Use Hungarian algorithm for graph matching
    cost_matrix = np.array(nx.to_numpy_array(graph_t2_reg))
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # Create temporal links based on matching results
    matching_reg_to_orig = [(i, j) for i, j in zip(range(len(center_2nd_reg)), col_ind) if cost_matrix[i, j] < np.inf]
    
    return matching_reg_to_orig

  

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

    #we first perform the registration of the images and the segmentations
    # path_warp, path_warpinv, path_second_image_reg, path_second_lesion_seg_reg = perform_registration(args.input_first_image, args.input_second_image, args.segmentation_first_image, 
    #                                                                                                     args.segmentation_second_image, args.output_folder)
    path_second_lesion_seg_reg = args.output_folder + '/second_lesion_seg_registered.nii.gz'
    
    #we then perform temporal lesion matching using the 1st time point lesion segmentation and the 2nd time point lesion segmentation registered to the 1st time point
    matching_table_1st_to_2nd_reg, center_1st, center_2nd_reg, first_data_labeled, second_data_reg_labeled = temporal_lesion_matching_on_reg_image(args.segmentation_first_image, path_second_lesion_seg_reg)
    
    #we then perfrom lesion clustering on the original lesion segmentation file from the second time point
    orig_second_lesion_seg = nib.load(args.segmentation_second_image)
    orig_second_lesion_seg_data = orig_second_lesion_seg.get_fdata()
    orig_second_lesion_seg_coordinates = np.argwhere(orig_second_lesion_seg_data == 1)
    clustering_2nd_orig = DBSCAN(eps=10, min_samples=5).fit(orig_second_lesion_seg_coordinates)
    orig_second_lesion_seg_labels = clustering_2nd_orig.labels_
    #we label the data by cluster
    orig_second_lesion_seg_data_labeled = np.zeros(orig_second_lesion_seg_data.shape)
    for i,voxel in enumerate(orig_second_lesion_seg_coordinates):
        orig_second_lesion_seg_data_labeled[voxel[0], voxel[1], voxel[2]] = orig_second_lesion_seg_labels[i]+1
    center_2nd_orig = {}
    #we compute the centers
    for lesion_2nd in np.unique(orig_second_lesion_seg_data_labeled):
        if lesion_2nd != 0:
            orig_second_lesion_seg_coordinates2 = np.argwhere(orig_second_lesion_seg_data_labeled == lesion_2nd)
            center_2nd_orig[int(lesion_2nd)] = np.mean(orig_second_lesion_seg_coordinates2, axis=0)


    """
    #we then perform lesion matching between the second time point lesion segmentation and the second time point lesion segmentation registered to the first time point
    second_reg_voxel_volume = np.prod(nib.load(args.input_first_image).header.get_zooms())
    second_voxel_volume = np.prod(nib.load(args.input_second_image).header.get_zooms())
    second_reg_volumes = volume_lesions(second_data_reg_labeled, second_reg_voxel_volume)
    second_orig_volumes = volume_lesions(orig_second_lesion_seg_data_labeled, second_voxel_volume)
    matching_reg_to_orig = graph_matching(center_2nd_orig, center_2nd_reg, second_orig_volumes, second_reg_volumes)


    
    #we build a matching table between the lesions of the first time point and the lesions of the second time point
    matching_table = np.empty(matching_table_1st_to_2nd_reg.shape)
    for i in range(len(matching_table_1st_to_2nd_reg)):
        matching_table[i][0] = matching_table_1st_to_2nd_reg[matching_reg_to_orig[i][1]][0]
        matching_table[i][1] = matching_table_1st_to_2nd_reg[matching_reg_to_orig[i][1]][1]
    
    #compute the volume of the lesions
    first_voxel_volume = np.prod(nib.load(args.input_first_image).header.get_zooms())
    second_voxel_volume = np.prod(nib.load(args.input_second_image).header.get_zooms())
    first_volumes = volume_lesions(first_data_labeled, first_voxel_volume)
    second_volumes = volume_lesions(orig_second_lesion_seg_data_labeled, second_voxel_volume)

    #we save the results in nifti files
    second_data_reg_labeled = nib.Nifti1Image(second_data_reg_labeled, orig_second_lesion_seg.affine, orig_second_lesion_seg.header)
    nib.save(second_data_reg_labeled, args.output_folder + '/second_lesion_seg_registered_clustered222.nii.gz')
    orig_second_lesion_seg_data_labeled = nib.Nifti1Image(orig_second_lesion_seg_data_labeled, orig_second_lesion_seg.affine, orig_second_lesion_seg.header)
    nib.save(orig_second_lesion_seg_data_labeled, args.output_folder + '/second_lesion_seg_clustered333.nii.gz')
    """

    #print output file
    """
    with open(args.output_folder + '/temporal_analysis_results.txt', 'w') as f:
        f.write('TEMPORAL ANALYSIS OF MULTIPLE SCLEROSIS LESIONS:\n')
        f.write('------------------------------------------------\n')
        f.write('Subject Name: at t0 ' + args.input_first_image.split('/')[-1] + ' and at t1 ' + args.input_second_image.split('/')[-1] + ' \n')
        f.write('------------------------------------------------\n')
        f.write('Lesions at t0: \n')
        f.write('------------------------------------------------\n')
        f.write('Number of lesions: ' + str(len(np.unique(first_data_labeled))-1) + '\n')
        for lesion in first_volumes:
            f.write('Lesion ' + str(int(lesion[0])) + ' volume: ' + str(lesion[1]) + ' mm3\n')
    """








if __name__ == '__main__':
    main()