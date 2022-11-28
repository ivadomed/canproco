# CanProCo

Code for preprocessing the CanProCo brain and spinal cord dataset

## Requirements

* SCT

## Dataset

Detailed instructions where data are stored and how to download them are available on the [intranet.neuro.polymtl.ca](https://intranet.neuro.polymtl.ca/computing-resources/data/git-datasets.html#usage).

## Preprocessing

1. Clone this repo

```commandline
git clone https://github.com/ivadomed/canproco.git
```

2. Run analysis across all subjects

```commandline
cd <PATH_TO_DATA>
sct_run_batch -c <PATH_TO_REPO>/etc/config_preprocess_data.json
```

Tip: You can run the analysis only across selected subjects. For details, see examples at the beginning of the preprocessing script.

Tip: Since analysis across many subjects can take a long time, it is recommended to run the analysis within a virtual terminal such as `screen`, details [here](https://intranet.neuro.polymtl.ca/geek-tips/bash-shell/README.html#screen-for-background-processes).

3. Manual correction

TODO

[Youtube manual](https://www.youtube.com/watch?v=lB-F8WOHGeg)