# CanProCo

Code for preprocessing the CanProCo brain and spinal cord dataset

## Requirements

* SCT
* Python 3
    * pyyaml
    * coloredlogs
* [ITK-SNAP](http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.SNAP3) for correcting cord segmentations

    **NOTE:** 
    Make sure to add ITK-SNAP to the system path:
    - For Windows, select the option during installation.
    - For macOS, after installation, go to **Help->Install Command-Line Tools**.

## Dataset

Detailed instructions where data are stored and how to download them are available on the [intranet.neuro.polymtl.ca](https://intranet.neuro.polymtl.ca/computing-resources/data/git-datasets.html#usage).
the dataset is under the name: `canproco`

## Preprocessing

### 1. Clone this repo

```commandline
git clone https://github.com/ivadomed/canproco.git
```

### 2. Run analysis across all subjects

```commandline
cd <PATH_TO_DATA>
sct_run_batch -c <PATH_TO_REPO>/etc/config_preprocess_data.json
```

Tip: You can run the analysis only across selected subjects. For details, see examples at the beginning of the preprocessing script.

Tip: Since analysis across many subjects can take a long time, it is recommended to run the analysis within a virtual terminal such as `screen`, details [here](https://intranet.neuro.polymtl.ca/geek-tips/bash-shell/README.html#screen-for-background-processes).

### 3. Manual correction of spinal cord segmentation

After running the analysis, check your Quality Control (qc) report by opening the file `./qc/index.html`. Use the "search" feature of the QC report to quickly jump to segmentations issues.

#### 1. Assess quality of segmentation

If segmentation issues are noticed while checking the quality report, proceed to manual correction using the procedure below:

1. In QC report, search for "deepseg" to only display results of spinal cord segmentation.
2. Review spinal cord segmentation.
3. Click on the `F` key to indicate if the segmentation/label is OK ✅, needs manual correction ❌ or if the data is not usable ⚠️ (artifact). Two .yml lists, one for manual corrections and one for the unusable data, will automatically be generated. 
4. Download the lists by clicking on `Download QC Fails` and on `Download Qc Artifacts`. 

The lists will have the following format:

*.yml list for correcting cord segmentation:*
~~~
FILES_SEG:
- sub-1000032_T1w.nii.gz
- sub-1000083_T2w.nii.gz
~~~

For the next steps, the script `/preprocess/manual_correction.py` loops through all the files listed in .yml file and opens an interactive window to manually correct segmentation. Each manually-corrected segmentation is saved under `derivatives/labels/` folder at the root of `PATH_DATA` according to the BIDS convention. Each manually-corrected file has the suffix `-manual`.

#### 2. Correct segmentations
For manual segmentation, you will need ITK-SNAP and this repository only. See **[Requirements](#requirements)**.

Here is a tutorial to manually correct segmentations. Note that the new QC report format with interactive features (✅/❌/⚠️) is not included in the tutorial.

[![IMAGE ALT TEXT](http://img.youtube.com/vi/vCVEGmKKY3o/sddefault.jpg)](https://youtu.be/vCVEGmKKY3o "Correcting segmentations across multiple subjects")

Run the following line and specify the .yml list for spinal cord segmentation with the flag `-config`:
~~~
cd canproco
python preprocess/manual_correction.py -config <.yml file> -path-in <PATH-PREPROCESSING>/data_processed -path-out <PATH_DATA>
~~~

After all corrections are done, you can generate a QC report by adding the flag `-qc-only-` to the command above. Note that SCT is required for generating QC report.

Here is another video for correcting manual segmentation: [Youtube manual](https://www.youtube.com/watch?v=lB-F8WOHGeg)

#### 3. Adding corrected segmentations to git-annex
After validating in the QC report the manual corrections, upload the manually-corrected segmentations in git-annex (see internal [documentation](https://intranet.neuro.polymtl.ca/computing-resources/data/git-datasets.html#upload)).

To update the dataset, add all manually-corrected files `derivatives/labels/`,  and include the qc zip file in the body of the PR.

## Manual lesion segmentation

Lesions are segmented on the PSIR or STIR contrasts (when available). The workflow is the following:

- Loop across subjects,
- Load PSIR/STIR in FSLeyes, and overlay the lesion segmentation (if it already exists),
- User create/adjust the segmentation mask using editing tools. Toggle the lesion overlay on/off to help validate the accuracy of the label ('Command F' on a Mac),
- Save the mask

We created a script to do this, which you can download [here](https://github.com/spinalcordtoolbox/manual-correction). 

Then, run:
```bash
python manual_correction.py -path-in <INPUT_PATH> -config <CONFIG_FILE> -path-out <OUTPUT_PATH> -viewer fsleyes
```

- `INPUT_PATH`: BIDS dataset from which manual lesions will be identified
- `CONFIG_FILE`: YML file that lists all subjects to be included in the manual lesion segmentation. This file can either be generated by SCT's QC report (see section "3. Manual correction of spinal cord segmentation") or manually. See example of a YML file below:
- `OUTPUT_PATH`: Optional. If not specified, output segmentations will be generated in the `INPUT_PATH`.

  <details><summary>Terminal output</summary>

  ```yaml
  FILES_LESION:
    - sub-edm005_ses-M0_PSIR.nii.gz
    - sub-edm008_ses-M0_PSIR.nii.gz
    - sub-edm010_ses-M0_PSIR.nii.gz
    - sub-edm011_ses-M0_PSIR.nii.gz
    - sub-edm013_ses-M0_PSIR.nii.gz
  ```

  </details>
