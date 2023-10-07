Segment SC using contrast-agnostic MONAI model from PSIR contrast and perform vertebral labeling

1. Create a conda virtual environment and install dependencies

```console
conda create -n monai python=3.8
conda activate monai
pip install -r sc_seg/requirements.txt
```

2. Segment SC using the contrast-agnostic MONAI model from PSIR contrast and perform vertebral labeling

The `sct_run_batch` wrapper script is used to run the `sc_seg/segment_sc.sh` across multiple subjects in parallel.

```console
sct_run_batch -config config.json
```

Example `config.json` file:

```json
{
  "path_data"   : "<PATH_TO_DATASET>/canproco",
  "path_output" : "<PATH_TO_DATASET>canproco_contrast-agnostic_2023-10-06",
  "script"      : "${HOME}/code/canproco/sc_seg/segment_sc.sh",
  "jobs"        : 1,
  "include"     : "ses-M0",
  "exclude"     : "sub-mon118",
  "script_args" : "${HOME}/code/canproco/sc_seg/run_inference_single_image.py ${HOME}/data/models/contrast-agnostic_final_monai_model/nnunet_nf=32_DS=1_opt=adam_lr=0.001_AdapW_CCrop_bs=2_64x192x320_20230918-2253"
}
```

Note that the `sc_seg/run_inference_single_image.py` script is just a copy of the 
[`contrast-agnostic-softseg-spinalcord/monai/run_inference_single_image.py`](https://github.com/sct-pipeline/contrast-agnostic-softseg-spinalcord/blob/nk/monai/monai/run_inference_single_image.py) script.

