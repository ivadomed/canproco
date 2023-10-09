The following shows the setup needed to run the `sc_seg_canproco.py`file`. 

- Create a conda environment 

```bash
conda create -n venv_monai python=3.9
```

- Activate the environment

```bash
conda activate venv_monai
```

- Install the requirements

```bash
pip install -r requirements_sc_seg.txt
```

## Running the spinal cord segmentation 

To run the spinal cord segmentation I used the following command: 

```bash
python sc_seg_canproco.py --input-data ~/canproco/project/data/canproco/  --contrast PSIR,STIR --qc-folder ~/canproco/project/canproco_final/qc_folder --seg-script ~/canproco/project/canproco_final/contrast-agnostic-softseg-spinalcord-nk-monai/monai/run_inference_single_image.py --path-to-model ~/duke/temp/muena/contrast-agnostic/final_monai_model/nnunet_nf=32_DS=1_opt=adam_lr=0.001_AdapW_CCrop_bs=2_64x192x320_20230918-2253 --output-folder ~/canproco/project/canproco_final/output_folder 
```


