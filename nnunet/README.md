# Installation instructions

## Installation of Anima metrics

Installation of 

```
cd ~
mkdir anima/
cd anima/
wget -q https://github.com/Inria-Empenn/Anima-Public/releases/download/v4.2/Anima-macOS-4.2.zip # for MACOS
unzip Anima-macOS-4.2.zip
rm Anima-macOS-4.2.zip
git lfs install
git clone --depth 1 https://github.com/Inria-Visages/Anima-Scripts-Public.git
git clone --depth 1 https://github.com/Inria-Visages/Anima-Scripts-Data-Public.git
```

Configure directories

```
cd ~
mkdir .anima/
touch .anima/config.txt

echo "[anima-scripts]" >> .anima/config.txt
echo "anima = ${HOME}/anima/Anima-Binaries-4.2/" >> .anima/config.txt
echo "anima-scripts-public-root = ${HOME}/anima/Anima-Scripts-Public/" >> .anima/config.txt
echo "extra-data-root = ${HOME}/anima/Anima-Scripts-Data-Public/" >> .anima/config.txt
```

## Installation of required libraries

Create a virtual invironment: 
~~~
conda create -n venv_nnunet python=3.9
~~~

Activate the environment with the following command:
~~~
conda activate venv_nnunet
~~~

To install required libraries to train an nnUNet v2:

```
pip install -r requirements_nnunet.txt
```

Install SpinalCordToolbox 6.0 :

Installation link : https://spinalcordtoolbox.com/user_section/installation.html


# Data preparation

Create the following folders:

~~~
mkdir nnUNet_raw
mkdir nnUNet_preprocessed
mkdir nnUNet_results
~~~

We are training a region-based nnUNet taking an image from contrasts PSIR or STIR and creating a mask with 0=background, 1=spinal cord and 2=MS lesion.

Convert the data to the nnUNet format :

~~~
python convert_BIDS_to_nnunet.py --path-data /path/to/BIDS/dataset --path-out /path/to/nnUNet_raw --taskname TASK-NAME --tasknumber DATASET-ID  --contrasts PSIR,STIR --test-ratio XX --time-point ses-XX --type training --exclude-file /path/to/exclude_file.yml
~~~

> **Note**
> The test ratio is 0.2 for 20% (train ratio is therefore 80%). For M0 images, the time point is ses-M0.

To mutliply PSIR images by -1 before training and convert the data to the nnUNet format :

~~~
python convert_BIDS_to_nnunet_with_mul_PSIR.py --path-data /path/to/BIDS/dataset --path-out /path/to/nnUNet_raw --taskname TASK-NAME --tasknumber DATASET-ID  --contrasts PSIR,STIR --test-ratio XX --time-point ses-XX --type training --exclude-file /path/to/exclude_file.yml
~~~

# Model training

Before training the model, nnU-Net performs data preprocessing and checks the integrity of the dataset:

~~~
export nnUNet_raw="/path/to/nnUNet_raw"
export nnUNet_preprocessed="/path/to/nnUNet_preprocessed"
export nnUNet_results="/path/to/nnUNet_results"

nnUNetv2_plan_and_preprocess -d DATASET-ID --verify_dataset_integrity
~~~

You will get the configuration plan for all four configurations (2d, 3d_fullres, 3d_lowres, 3d_cascade_fullres).

To train the model, use the following command:
~~~
CUDA_VISIBLE_DEVICES=XXX nnUNetv2_train DATASET-ID CONFIG FOLD --npz
~~~

> **Note**
> Example for Dataset 101, on 2d config on fold 0: CUDA_VISIBLE_DEVICES=2 nnUNetv2_train 101 2d 0 --npz

# Model inference

Convert data to nnUNet format for inference using `convert_BIDS_to_nnunet.py` or `convert_BIDS_to_nnunet_with_mul_PSIR.py` with `--type=inference`.

Then perform inference:
~~~
CUDA_VISIBLE_DEVICES=XXX nnUNetv2_predict -i /path/to/image/folder -o /path/to/predictions -d DATASET_ID -c CONFIG --save_probabilities -chk checkpoint_best.pth -f FOLD
~~~

# Inference evaluation

First, convert the predictions back the BIDS format ; this only keeps the lesion segmentation and discards spinal cord segmentation :

~~~
python convert_predictions_to_BIDS.py --pred-folder /path/to/predictions --out-folder /path/to/output/folder --conversion-dict /path/to/conversion/dict
~~~

If you are converting predictions which are not the evaluations set from nnUNet, use the flag `--not-imageTs`.

Then, you can evaluate the lesion prediction with Anima metrics

~~~
python evaluate_lesion_seg_prediction.py --pred-folder path/to/predictions --dataset path/to/dataset --animaPath path/to/animaSegPerfAnalyzer --output-folder path/to/output_folder
~~~

# Evaluation analysis

The following Notebook `nnUNet_inference_analysis.ipynb` was used to perform analysis of the nnUNet segmentations. 
To use it, change the path to the csv files. 
