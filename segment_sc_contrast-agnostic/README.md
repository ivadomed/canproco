# Segment SC using contrast-agnostic MONAI model from STIR/PSIR contrast and perform vertebral labeling

This part of the repo was created to segment the spinal cord of every subject in CanProCo. Below, we detail a recommended method (latest) and the one used previously. 

## New method to segment spinal cord (recommended)

This method uses the Spinal Cord Toolbox (SCT)

### Install dependencies

- [Spinal Cord Toolbox (SCT) v6.2](https://github.com/spinalcordtoolbox/spinalcordtoolbox/releases/tag/6.2) or higher -- follow the installation instructions [here](https://github.com/spinalcordtoolbox/spinalcordtoolbox?tab=readme-ov-file#installation)
- [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) 
- Python

Once the dependencies are installed, download the latest rootlets model:

```bash
sct_deepseg -install-task seg_sc_contrast_agnostic
```

### Getting the rootlet segmentation

To segment a single image, run the following command: 

```bash
sct_deepseg -i <INPUT> -o <OUTPUT> -task seg_sc_contrast_agnostic
```

For example:

```bash
sct_deepseg -i sub-001_T2w.nii.gz -o sub-001_T2w_sc-seg.nii.gz -task seg_sc_contrast_agnostic
```

## Our previous method of segmenting the spinal cord (not recommended)

1. Create a conda virtual environment and install dependencies

```console
yes | conda create -n monai python=3.9
conda activate monai
yes | pip install -r segment_sc_contrast-agnostic/requirements.txt
```

2. Segment SC using the contrast-agnostic MONAI model from STIR/PSIR contrast and perform vertebral labeling

The `sct_run_batch` wrapper script is used to run the `segment_sc_contrast-agnostic/segment_sc_contrast-agnostic.sh` across multiple subjects in parallel.

```console
sct_run_batch -config config.json
```

Example `config.json` file:

```json
{
  "path_data"   : "<PATH_TO_DATASET>/canproco",
  "path_output" : "<PATH_TO_DATASET>canproco_contrast-agnostic_2023-10-06",
  "script"      : "${HOME}/code/canproco/segment_sc_contrast-agnostic/01_segment_sc_contrast-agnostic.sh",
  "jobs"        : 1,
  "include"     : "ses-M0",
  "exclude_list": "sub-mon118 sub-mon006 sub-mon009 sub-mon032 sub-mon097 sub-mon148 sub-mon168 sub-mon191 sub-van176 sub-van206 sub-tor133 sub-cal149",
  "script_args" : "${HOME}/code/canproco/segment_sc_contrast-agnostic/run_inference_single_image.py ${HOME}/data/models/contrast-agnostic_final_monai_model/nnunet_nf=32_DS=1_opt=adam_lr=0.001_AdapW_CCrop_bs=2_64x192x320_20230918-2253"
}
```

ℹ️ The `exclude_list` is a list of subjects to exclude from the batch processing. It is based on the https://github.com/ivadomed/canproco/blob/main/exclude.yml file.

ℹ️ The conda environment with MONAI is required to run the `segment_sc_contrast-agnostic/segment_sc_contrast-agnostic.sh` script.

ℹ️ Note that the `segment_sc_contrast-agnostic/run_inference_single_image.py` script is just a copy of the 
[`contrast-agnostic-softseg-spinalcord/monai/run_inference_single_image.py`](https://github.com/sct-pipeline/contrast-agnostic-softseg-spinalcord/blob/nk/monai/monai/run_inference_single_image.py) script.
