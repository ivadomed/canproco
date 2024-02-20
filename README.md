# CanProCo

This repository contains the code used in the different projects of the Canadian Prospective Cohort to Understand Progression in People Living with MS (CanProCo). 


For manual correction and manual segmentation of MS lesions, please refer to this repository : [manual-correction](https://github.com/spinalcordtoolbox/manual-correction).


## Table of Contents
=================
* [Computing Cross-Sectional Area (CSA)](https://github.com/ivadomed/canproco/blob/main/scripts-t2w_csa/README.md)
* [Spinal-cord segmentation](https://github.com/ivadomed/canproco/blob/main/segment_sc_contrast-agnostic/README.md)
* [nnUNet training](https://github.com/ivadomed/canproco/blob/main/nnunet/README.md)
* [Dataset analysis of MS lesions](https://github.com/ivadomed/canproco/blob/main/dataset_analysis/README.md)
* [Building lesion frequency maps](https://github.com/ivadomed/canproco/blob/main/lesion-mapping/README.md)
* [Using the model for MS lesion segmentation](#how-to-use-the-model-for-ms-lesion-segmentation)

## How to use the model for MS lesion segmentation

### Step 1: Cloning the Repository

Open a terminal and clone the repository using the following command:

~~~
git clone https://github.com/ivadomed/canproco.git
~~~

### Step 2: Setting up the Environment

The following commands show how to set up the environment. Note that the documentation assumes that the user has `conda` installed on their system. Instructions on installing `conda` can be found [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

1. Create a conda environment with the following command:
```
conda create -n venv_nnunet python=3.9
```

2. Activate the environment with the following command:
```
conda activate venv_nnunet
```

3. Install the required packages with the following command:
```
cd canproco
pip install -r packaging/requirements.txt
```

4. Download the model (`model_ms_seg_sc-lesion_regionBased.zip`) from the repository's latest release, which can be found [here](https://github.com/ivadomed/canproco/releases/tag/r20240125), and unzip it. 
 
### Step 3: Getting the Predictions

To segment a single image using the trained model, run the following command from the terminal. This assumes that the model has been downloaded and is available locally. The release contains two models, the 2D nnUNet as well as the 3D nnUNet. Our experiments showed that both worked similarly. 

```bash
python packaging/run_inference_single_subject.py --path-image /path/to/image --path-out /path/to/output/directory --path-model /path/to/model 
```

The output contains the spinal cord segmentation (with value 1) and the MS lesion segmentation (with value 2). It uses a region-based approach, meaning that lesions are always located within the spinal cord segmentation. 

ℹ️ The script also supports getting segmentations on a GPU. To do so, simply add the flag `--use-gpu` at the end of the above commands. By default, the inference is run on the CPU. It is useful to note that obtaining the predictions from the GPU is significantly faster than the CPU.
