# Lesion mapping

### 1. Register STIR/PSIR images to the PAM50 space 

Register STIR/PSIR images to the PAM50 space and bring the GT lesion and spinal cord masks (located under 
`derivatives/labels`) to the PAM50 template space.

```console
sct_run_batch -config config.json
```

Example `config.json` file:

```json
{
 "path_data"   : "<PATH_TO_DATASET>/canproco",
 "path_output" : "<PATH_TO_DATASET>/canproco_register_to_PAM50_2023-10-21",
 "script"      : "<PATH_TO_REPO>/canproco/lesion-mapping/01_register_to_pam50.sh",
 "jobs"        : 16,
 "exclude_list": "sub-mon118 sub-mon006 ..." 
}
```

### 2. Generate the lesion frequency map

Generate the lesion frequency map (LFM) in the PAM50 space.

ℹ️ The `02_generate_lesion_frequency_maps.py` script requires the SCT conda environment to be activated:
```console
source ${SCT_DIR}/python/etc/profile.d/conda.sh
conda activate venv_sct
```

Run the script:

```console
python 02_generate_lesion_frequency_maps.py \
    -ifolder canproco_register_to_PAM50_2023-10-21/data_processed \
    -participants-tsv canproco/participants.tsv \
    -exclude-yml exclude_M0_M12_comparison.yml \
    -ofolder canproco_register_to_PAM50_2023-10-21/results \
    -session M0
```

### 3. Generate png images of the LFM

From the Lesion Frequency Map (LFM), generate:
    - axial png images for each vertebral level (average across axial slices for each vertebral level)
    - a single sagittal png image (average across sagittal slices)

ℹ️ The `03_lfm_to_png.py` script requires the SCT conda environment to be activated:
```console
source ${SCT_DIR}/python/etc/profile.d/conda.sh
conda activate venv_sct
```

Run the script:

```console
python 03_lfm_to_png.py \
    -lfm-path spinalcord_LFM_MS_200_participants.nii.gz \
    -ofolder lfm \
    -thr 0.15 \
    -include-gm
```