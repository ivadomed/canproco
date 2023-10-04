"""
This python file is used to analyse the evolution of the lesions over time.
It first performs a registration of the images and their segmentation on the two time points.
Then, it performs clustering of the lesions on each image and computes the volume and the center of each lesion. 
Finally, it compares the lesions between the two time points and computes the evolution of the lesions.
This file relies on spinal cord toolbox: version 6.0

Args:
    -i-1, --input-first-image: path to the first image 
    -i-2, --input-second-image: path to the second image
    -seg-1, --segmentation-first-image: path to the lesion segmentation of the first image
    -seg-2, --segmentation-second-image: path to the lesion segmentation of the second image
    -sc-seg-1, --sc-segmentation-first-image: path to the spinal cord segmentation of the first image
    -sc-seg-2, --sc-segmentation-second-image: path to the spinal cord segmentation of the second image
    -o, --output-folder: path to the output folder

Returns:
    None

Example:
    python time_point_lesion_evolution.py -i-1 path/to/first/image -i-2 path/to/second/image -seg-1 path/to/first/segmentation -seg-2 path/to/second/segmentation -o path/to/output/folder

To do:
    * Remove some of the outputed files

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
from time import time


def get_parser():
    """
    This function parses the arguments given to the script.

    Args:
        None

    Returns:
        parser: parser containing the arguments
    """

    parser = argparse.ArgumentParser(description='Analyse and compare the results of the segmentation of the MS lesion on the spinal cord at two time point.')
    parser.add_argument('-i-1', '--input-first-image', required=True,
                        help='Path to the first image')
    parser.add_argument('-i-2', '--input-second-image', required=True,
                        help='Path to the second image')
    parser.add_argument('-seg-1', '--segmentation-first-image', required=True,
                        help='Path to the lesion segmentation of the first image')
    parser.add_argument('-seg-2', '--segmentation-second-image', required=True,
                        help='Path to the lesion segmentation of the second image')
    parser.add_argument ('-sc-seg-1', '--sc-segmentation-first-image', required=False, default=None,
                        help='Path to the spinal cord segmentation of the first image')
    parser.add_argument ('-sc-seg-2', '--sc-segmentation-second-image', required=False, default=None,
                        help='Path to the spinal cord segmentation of the second image')
    parser.add_argument('-o', '--output-folder', required=True,
                        help='Path to the output folder')
    
    return parser


def perform_registration(path_first_image, path_second_image, path_first_lesion_seg, path_second_lesion_seg, path_sc_seg_first_image, path_sc_seg_second_image, path_output_folder):
    """
    This function performs the registration of the images and the segmentations using spinal cord toolbox

    Args:
        path_first_image: path to the first image
        path_second_image: path to the second image
        path_first_lesion_seg: path to the lesion segmentation of the first image
        path_second_lesion_seg: path to the lesion segmentation of the second image
        path_sc_seg_first_image: path to the spinal cord segmentation of the first image
        path_sc_seg_second_image: path to the spinal cord segmentation of the second image
        path_output_folder: path to the output folder
    
    Returns:
        path_warp: path to the warp file
        path_warpinv: path to the inverse warp file
        path_second_image_registered: path to the second image registered
        path_second_lesion_seg_registered: path to the second lesion segmentation registered
    """

    # Segmentation of the spinal cord on the two images
    ## If the spinal cord segmentation is given, we use it:
    if path_sc_seg_first_image is not None:
        first_sc_output_path = path_sc_seg_first_image
    ## Otherwise, we compute it using sct_deepseg_sc
    else:
        first_sc_output_path = os.path.join(path_output_folder, 'first_image_sc_seg.nii.gz')
        os.system('sct_deepseg_sc -i ' + path_first_image + ' -c t2 -o ' + first_sc_output_path + ' -qc ' + os.path.join(path_output_folder, 'QC_folder'))

    ## Same for the second image
    if path_sc_seg_second_image is not None:
        second_sc_output_path = path_sc_seg_second_image
    else:
        second_sc_output_path = os.path.join(path_output_folder, 'second_image_sc_seg.nii.gz')
        os.system('sct_deepseg_sc -i ' + path_second_image + ' -c t2 -o ' + second_sc_output_path + ' -qc ' + os.path.join(path_output_folder, 'QC_folder'))    

    # Then we compute the vertebral levels of the images using sct_label_vertebrae
    first_vert_seg_path = os.path.join(path_output_folder, 'first_image_sc_seg.nii.gz')
    second_vert_seg_path = os.path.join(path_output_folder, 'second_image_sc_seg.nii.gz')
    os.system('sct_label_vertebrae -i ' + path_first_image + ' -s ' + first_vert_seg_path + ' -c t2 -ofolder ' + path_output_folder + ' -qc ' + os.path.join(path_output_folder, 'QC_folder'))
    os.system('sct_label_vertebrae -i ' + path_second_image + ' -s ' + second_vert_seg_path + ' -c t2 -ofolder ' + path_output_folder + ' -qc ' + os.path.join(path_output_folder, 'QC_folder'))

    # We remove the vertebral levels that are not in common in both images using 'sct_label_utils -remove-reference'
    first_vert_output_path = os.path.join(path_output_folder, 'first_image_sc_seg_labeled.nii.gz')
    second_vert_output_path = os.path.join(path_output_folder, 'second_image_sc_seg_labeled.nii.gz')
    t0=time()
    print(t0)
    os.system('sct_label_utils -i ' + first_vert_output_path + ' -remove-reference ' + second_vert_output_path + ' -o ' + first_vert_output_path)
    t1 = time()
    print("time to remove-reference",t1-t0)
    os.system('sct_label_utils -i ' + second_vert_output_path + ' -remove-reference ' + first_vert_output_path + ' -o ' + second_vert_output_path)
    t2 = time()
    print("time to remove-reference",t2-t1)

    """
    ## We first get the labels of the two images
    first_label = nib.load(os.path.join(path_output_folder, 'first_image_sc_seg_labeled.nii.gz'))
    second_label = nib.load(os.path.join(path_output_folder, 'second_image_sc_seg_labeled.nii.gz'))
    first_label_data = first_label.get_fdata()
    second_label_data = second_label.get_fdata()
    common_labels = np.intersect1d(np.unique(first_label_data), np.unique(second_label_data))

    # We then remove the labels that are not in common
    for label in np.unique(first_label_data):
        if label not in common_labels:
            first_label_data[first_label_data == label] = 0
    # We overwrite the segmentation to keep only common labels
    first_label = nib.Nifti1Image(first_label_data, first_label.affine, first_label.header)
    nib.save(first_label, os.path.join(path_output_folder, 'first_image_sc_seg_labeled.nii.gz'))
    # Same for the second image
    for label in np.unique(second_label_data):
        if label not in common_labels:
            second_label_data[second_label_data == label] = 0
    second_label = nib.Nifti1Image(second_label_data, second_label.affine, second_label.header)
    nib.save(second_label, os.path.join(path_output_folder, 'second_image_sc_seg_labeled.nii.gz'))
    
    """

    # Then we register the images and the segmentations with the first image as reference using sct_register_multimodal using both the spinal cord segmentation and the vertebral levels
    parameters = 'step=0,type=label,dof=Tx_Ty_Tz:step=1,type=im,algo=dl'
    os.system('sct_register_multimodal -i ' + path_second_image + ' -d ' + path_first_image + ' -ilabel ' + os.path.join(path_output_folder,'second_image_sc_seg_labeled.nii.gz') 
                + ' -dlabel ' + os.path.join(path_output_folder,'first_image_sc_seg_labeled.nii.gz') + ' -dseg ' + first_sc_output_path + ' -ofolder ' + path_output_folder + ' -owarp ' + os.path.join(path_output_folder,'warp_2nd_to_1st.nii.gz')
                 + ' -owarpinv ' + os.path.join(path_output_folder, 'warp_1st_to_2nd.nii.gz') + ' -param ' + parameters + ' -x linear -qc ' + os.path.join(path_output_folder, 'QC_folder'))

    # We apply the warping to the segmentation of the second image
    os.system('sct_apply_transfo -i ' + path_second_lesion_seg + ' -d ' + path_first_lesion_seg + ' -w ' + os.path.join(path_output_folder,'warp_2nd_to_1st.nii.gz') 
              + ' -o ' + os.path.join(path_output_folder, 'second_lesion_seg_registered.nii.gz') + ' -x nn')

    return os.path.join(path_output_folder,'warp_2nd_to_1st.nii.gz'), os.path.join(path_output_folder, 'warp_1st_to_2nd.nii.gz'), os.path.join(path_output_folder, 'second_image_sc_seg_registered.nii.gz'), os.path.join(path_output_folder, 'second_lesion_seg_registered.nii.gz')


def temporal_lesion_matching_on_reg_image(path_first_lesion_seg, path_second_lesion_seg_reg, path_output_folder):
    """
    This function performs temporal lesion matching using the 1st time point lesion segmentation and the 2nd time point lesion segmentation registered to the 1st time point.
    The matching is done using overlapping. If more than 50% of the lesion of the 1st time point is overlapping with a lesion of the 2nd time point,
    then we consider that the two lesions are the same.

    Args:
        path_first_lesion_seg: path to the lesion segmentation of the first image
        path_second_lesion_seg_reg: path to the lesion segmentation of the second image registered to the first image
        path_output_folder: path to the output folder

    Returns:
        matching_table: table containing the matching between lesions of the first time point and lesions of the second time point
        center_1st: coordinates of the center of the lesions of the first time point
        center_2nd: coordinates of the center of the lesions of the second time point
        1st_data_labeled: segmentation of the first time point with each lesion in a different color
        2nd_data_labeled: segmentation of the second time point with each lesion in a different color
    """

    # Just in case, we binarize the second lesion segmentation registered to the first image to avoid any problem
    second_lesion_seg_reg = nib.load(path_second_lesion_seg_reg)
    second_lesion_seg_reg_data = second_lesion_seg_reg.get_fdata()
    second_lesion_seg_reg_data[second_lesion_seg_reg_data > 1e-5] = 1

    # We then load data from the first lesion segmentation
    first_lesion_seg = nib.load(path_first_lesion_seg)
    first_lesion_seg_data = first_lesion_seg.get_fdata()

    # We then perform clustering on the two images
    # We first get the coordinates of the lesions
    first_lesion_seg_coordinates = np.argwhere(first_lesion_seg_data == 1)
    second_lesion_seg_reg_coordinates = np.argwhere(second_lesion_seg_reg_data == 1)
    # We then perform clustering on the coordinates
    clustering_1st = DBSCAN(eps=10, min_samples=5).fit(first_lesion_seg_coordinates)
    first_lesion_seg_labels = clustering_1st.labels_
    clustering_2nd = DBSCAN(eps=10, min_samples=5).fit(second_lesion_seg_reg_coordinates)
    second_lesion_seg_reg_labels = clustering_2nd.labels_

    # Create the same tables with the lesions with unique voxel value (voxel value = lesion ID)
    first_lesion_seg_data2 = np.zeros(first_lesion_seg_data.shape)
    second_lesion_seg_reg_data2 = np.zeros(second_lesion_seg_reg_data.shape)
    for i,voxel in enumerate(first_lesion_seg_coordinates):
        # Adding "1" to first_lesion_seg_labels[i] because the first label is 0.
        first_lesion_seg_data2[voxel[0], voxel[1], voxel[2]] = first_lesion_seg_labels[i]+1
    for i,voxel in enumerate(second_lesion_seg_reg_coordinates):
        # Same here
        second_lesion_seg_reg_data2[voxel[0], voxel[1], voxel[2]] = second_lesion_seg_reg_labels[i]+1

    # For each lesion of the second time point we consider that it corresponds to the lesion of the first time point with which it overlaps more than 50%
    list_lesion_2nd = []
    corresponding_lesion_1st = []
    overlapping_ratio = []
    center_2nd = []
    center_1st = []
    first_pass = True
    for lesion_2nd in np.unique(second_lesion_seg_reg_data2):
        if lesion_2nd != 0:
            list_lesion_2nd.append(lesion_2nd)
            # We compute the centers
            second_lesion_iter_center = np.argwhere(second_lesion_seg_reg_data2 == lesion_2nd)
            center_2nd.append(np.mean(second_lesion_iter_center, axis=0))
            for lesion_1st in np.unique(first_lesion_seg_data2):
                # We compute the centers
                if first_pass:
                    first_lesion_iter_center = np.argwhere(first_lesion_seg_data2 == lesion_1st)
                    center_1st.append(np.mean(first_lesion_iter_center, axis=0))
                    first_pass = False
                # We compute the overlap
                if lesion_1st !=0:
                    # We get the coordinates of the lesion
                    second_lesion_seg_coordinates2 = np.argwhere(second_lesion_seg_reg_data2 == lesion_2nd)
                    # We get the coordinates of the lesion in the first time point
                    first_lesion_seg_coordinates = np.argwhere(first_lesion_seg_data2 == lesion_1st)
                    # We compute the overlap of the two sets of coordinates
                    second_set = set( tuple(points) for points in second_lesion_seg_coordinates2)
                    first_set = set(tuple(points) for points in  first_lesion_seg_coordinates)
                    intersection = len(second_set.intersection(first_set))
                    ratio = intersection/len(first_set)
                    if ratio > 0.5:
                        corresponding_lesion_1st.append(lesion_1st)
                        overlapping_ratio.append(ratio)
                    else:
                        corresponding_lesion_1st.append(None)
                        overlapping_ratio.append(None)

    lesion_dict = {"list_lesion_2nd": list_lesion_2nd, "corresponding_lesion_1st": corresponding_lesion_1st, "overlapping_ratio": overlapping_ratio,
                    "center_2nd": center_2nd, "center_1st": center_1st}

    # Finally we overwrite the segmentation files with the lesions in colors matching from t1 to t2
    final_first_lesion_seg = np.zeros(first_lesion_seg_data2.shape)
    for lesion in np.unique(first_lesion_seg_data2):
        if lesion!=0:
            arg_where_lesion = np.argwhere(first_lesion_seg_data2 == lesion)
            index_lesion = np.argwhere(lesion_dict['corresponding_lesion_1st'] == lesion)[0]
            for coord in arg_where_lesion:
                final_first_lesion_seg[coord[0], coord[1], coord[2]] = lesion_dict['list_lesion_2nd'][index_lesion[0]]
    final_second_lesion_seg_reg = np.copy(second_lesion_seg_reg_data2)

    # We save the results in nifti files
    first_lesion_seg2 = nib.Nifti1Image(first_lesion_seg_data2, first_lesion_seg.affine, first_lesion_seg.header)
    second_lesion_seg2 = nib.Nifti1Image(final_second_lesion_seg_reg, second_lesion_seg_reg.affine, second_lesion_seg_reg.header)
    nib.save(first_lesion_seg2, os.path.join(path_output_folder, 'first_lesion_seg_clustered.nii.gz'))
    nib.save(second_lesion_seg2, os.path.join(path_output_folder, 'second_lesion_seg_clustered.nii.gz'))
    
    return lesion_dict, final_first_lesion_seg, final_second_lesion_seg_reg


def lesion_matching(path_t2_seg, path_t2_inv_warpped, path_output_folder):
    """
    This function performs lesion matching between the 2nd time point lesion segmentation and the 2nd time point lesion segmentation inversely warpped.

    Args:
        path_t2_seg: path to the lesion segmentation of the second time point
        path_t2_inv_warpped: path to the lesion segmentation of the second time point inversely warpped
        path_output_folder: path to the output folder

    Returns:
        path_t2_seg: path to the lesion segmentation of the second time point with the lesions in same color as the lesions of the second time point inversely warpped
    """
    # We first load data
    t2_seg = nib.load(path_t2_seg)
    t2_seg_data = t2_seg.get_fdata()
    t2_inv_warpped = nib.load(path_t2_inv_warpped)
    t2_inv_warpped_data = t2_inv_warpped.get_fdata()

    # We then perform clustering on t2_seg
    t2_seg_coordinates = np.argwhere(t2_seg_data == 1)
    clustering_t2 = DBSCAN(eps=10, min_samples=5).fit(t2_seg_coordinates)
    t2_seg_labels = clustering_t2.labels_
    # From this we build a dataset with the lesions in different colors
    t2_seg_data_labeled = np.zeros(t2_seg_data.shape)
    for i,voxel in enumerate(t2_seg_coordinates):
        t2_seg_data_labeled[voxel[0], voxel[1], voxel[2]] = t2_seg_labels[i]+1
    
    # We now compute the centers of the lesions
    center_t2 = []
    color_t2 = []
    for lesion_t2 in np.unique(t2_seg_data_labeled):
        if lesion_t2 != 0:
            t2_seg_coordinates2 = np.argwhere(t2_seg_data_labeled == lesion_t2)
            center_t2.append(np.mean(t2_seg_coordinates2, axis=0))
            color_t2.append(lesion_t2)
    center_t2_inv_warpped = []
    color_t2_inv_warpped = []
    for lesion_t2 in np.unique(t2_inv_warpped_data):
        if lesion_t2 != 0:
            t2_inv_warpped_coordinates = np.argwhere(t2_inv_warpped_data == lesion_t2)
            center_t2_inv_warpped.append(np.mean(t2_inv_warpped_coordinates, axis=0))
            color_t2_inv_warpped.append(lesion_t2)

    if len(center_t2) != len(center_t2_inv_warpped):
        print('Error: the number of lesions is different between the two segmentations.')
        return None

    # We then perform matching by closest distance of the two center lists
    matching_table = np.zeros((len(center_t2), 2))
    for i, center in enumerate(center_t2):
        min_distance = np.inf
        for j, center_inv_warpped in enumerate(center_t2_inv_warpped):
            distance = np.sqrt((center[0] - center_inv_warpped[0])**2 + (center[1] - center_inv_warpped[1])**2 + (center[2] - center_inv_warpped[2])**2)
            if distance < min_distance:
                min_distance = distance
                matching_table[i][0] = i
                matching_table[i][1] = j
    
    # We then overwrite the segmentation file with the lesions in colors matching from t1 to t2
    final_t2_seg = np.zeros(t2_seg_data_labeled.shape)
    for i,lesion_value in enumerate(color_t2):
        final_t2_seg[t2_seg_data_labeled == lesion_value] = color_t2_inv_warpped[int(matching_table[i][1])]
    final_t2_seg = nib.Nifti1Image(final_t2_seg, t2_seg.affine, t2_seg.header)
    out_path = path_output_folder + path_t2_seg.split('/')[-1].split('.')[0] + '_matched.nii.gz'
    nib.save(final_t2_seg, out_path)

    return out_path
    

def comparison_lesions_t1_t2(path_first_lesion_seg, path_second_lesion_seg):
    """
    This function performs comparison of lesions from t1 and t2. 

    Args:
        path_first_lesion_seg: path to the lesion segmentation of the first image
        path_second_lesion_seg: path to the lesion segmentation of the second image

    Returns:
        t1_dict: dictionary containing the list of lesions, the volume of the lesions and the center of the lesions for the first time point
        t2_dict: dictionary containing the list of lesions, the volume of the lesions and the center of the lesions for the second time point
    """

    # We first load data
    first_lesion_seg = nib.load(path_first_lesion_seg)
    first_lesion_seg_data = first_lesion_seg.get_fdata()
    second_lesion_seg = nib.load(path_second_lesion_seg)
    second_lesion_seg_data = second_lesion_seg.get_fdata()

    # We get voxel volume
    voxel_volume_t1 = np.prod(first_lesion_seg.header.get_zooms())
    voxel_volume_t2 = np.prod(second_lesion_seg.header.get_zooms())

    # We iterate over the lesions of the first time point
    list_lesions_t1 = np.unique(first_lesion_seg_data)
    list_lesions_t1 = list_lesions_t1[1:] #removed 0
    volume_lesions_t1 = []
    center_t1 = []
    for lesion in list_lesions_t1:
        volume_lesions_t1.append(len(np.argwhere(first_lesion_seg_data == lesion))*voxel_volume_t1)
        center_t1.append(np.mean(np.argwhere(first_lesion_seg_data == lesion), axis=0))
    
    rounded_list_t1 = [round(elem) for elem in list_lesions_t1]
    t1_dict = {"list_lesions": rounded_list_t1, "volume_lesions": volume_lesions_t1, "center": center_t1}
    
    # We iterate over the lesions of the second time point
    list_lesions_t2 = np.unique(second_lesion_seg_data)
    list_lesions_t2 = list_lesions_t2[1:] #removed 0
    volume_lesions_t2 = []
    center_t2 = []
    for lesion in list_lesions_t2:
        volume_lesions_t2.append(len(np.argwhere(second_lesion_seg_data == lesion))*voxel_volume_t2)
        center_t2.append(np.mean(np.argwhere(second_lesion_seg_data == lesion), axis=0))
    
    rounded_list_t2 = [round(elem) for elem in list_lesions_t2]
    t2_dict = {"list_lesions": rounded_list_t2, "volume_lesions": volume_lesions_t2, "center": center_t2}

    return t1_dict, t2_dict
    

def main():
    """
    This function is the main function of the script.
    
    Args:
        None

    Returns:
        None
    """

    # Get the parser
    parser = get_parser()
    args = parser.parse_args()

    # We first perform the registration of the images and the segmentations
    path_warp, path_warpinv, path_second_image_reg, path_second_lesion_seg_reg = perform_registration(args.input_first_image, args.input_second_image, args.segmentation_first_image, 
                                                                                                         args.segmentation_second_image, args.sc_segmentation_first_image,
                                                                                                         args.sc_segmentation_second_image, args.output_folder)
    
    # We then perform temporal lesion matching using the 1st time point lesion segmentation and the 2nd time point lesion segmentation registered to the 1st time point
    lesion_dict, final_first_lesion_seg, final_second_lesion_seg_reg = temporal_lesion_matching_on_reg_image(args.segmentation_first_image, path_second_lesion_seg_reg, args.output_folder)
    
    # We then perform inverse warpping on the image 
    os.system('sct_apply_transfo -i ' + os.path.join(args.output_folder, 'second_lesion_seg_clustered.nii.gz') + ' -d ' + args.segmentation_second_image + ' -w ' + os.path.join(args.output_folder, 'warp_1st_to_2nd.nii.gz') 
              + ' -o ' + os.path.join(args.output_folder, 'second_lesion_seg_inv_warpped.nii.gz') + ' -x nn')
    
    # Now we perform lesion matching between the 2nd time point lesion segmentation and the 2nd time point lesion segmentation inversely warpped
    t2_final_seg = lesion_matching(args.segmentation_second_image, os.path.join(args.output_folder, 'second_lesion_seg_inv_warpped.nii.gz'), args.output_folder)

    # Now we perform comparison of lesions from t1 and t2
    t1_dict, t2_dict = comparison_lesions_t1_t2(args.segmentation_first_image,t2_final_seg)

    # We then print the output file
    list_new_lesions = set(t2_dict['list_lesions']).difference(set(t1_dict['list_lesions']))
    common_lesions = list(set(t2_dict['list_lesions']).intersection(set(t1_dict['list_lesions'])))
    
    with open(args.output_folder + '/temporal_analysis_results.txt', 'w') as f:
        f.write('TEMPORAL ANALYSIS OF MULTIPLE SCLEROSIS LESIONS:\n')
        f.write('------------------------------------------------\n')
        f.write('Subject Name: at ses-0 ' + args.input_first_image.split('/')[-1] + ' and at ses-1 ' + args.input_second_image.split('/')[-1] + ' \n')
        f.write('------------------------------------------------\n')
        f.write('Lesions at ses-0: ' + str(len(t1_dict['list_lesions'])) + '\n')
        f.write('Lesions at ses-1: ' + str(len(t2_dict['list_lesions'])) + '\n')
        f.write('------------------------------------------------\n')
        f.write('Lesions at ses-0: \n')
        for lesion in t1_dict['list_lesions']:
            f.write('Lesion ' + str(int(lesion)) + ' \n')
            f.write('Volume: ' + str(t1_dict['volume_lesions'][t1_dict['list_lesions'].index(lesion)]) + ' mm3\n')
            f.write('Center: ' + str(t1_dict['center'][t1_dict['list_lesions'].index(lesion)]) + '\n')
            f.write('\n')
        f.write('------------------------------------------------\n')
        f.write('Lesions at ses-1: \n')
        for lesion in t2_dict['list_lesions']:
            f.write('Lesion ' + str(int(lesion)) + ' \n')
            f.write('Volume: ' + str(t2_dict['volume_lesions'][t2_dict['list_lesions'].index(lesion)]) + ' mm3\n')
            f.write('Center: ' + str(t2_dict['center'][t2_dict['list_lesions'].index(lesion)]) + '\n')
            f.write('\n')
        f.write('------------------------------------------------\n')
        f.write('EVOLUTION:')
        for lesion in common_lesions:
            f.write('Lesion ' + str(int(lesion)) + ' \n')
            f.write('Volume at ses-0: ' + str(t1_dict['volume_lesions'][t1_dict['list_lesions'].index(lesion)]) + ' mm3\n')
            f.write('Volume at ses-1: ' + str(t2_dict['volume_lesions'][t2_dict['list_lesions'].index(lesion)]) + ' mm3\n')
            f.write('Volume increase in %: ' + str((t2_dict['volume_lesions'][t2_dict['list_lesions'].index(lesion)] - t1_dict['volume_lesions'][t1_dict['list_lesions'].index(lesion)])/t1_dict['volume_lesions'][t1_dict['list_lesions'].index(lesion)]*100) + ' %\n')
        f.write('------------------------------------------------\n')


if __name__ == '__main__':
    main()