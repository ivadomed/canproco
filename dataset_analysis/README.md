# Installation

To install the required libraries run :

~~~
pip install -r requirements_dataframe.txt
~~~

It is also required to install SpinalCordToolbox 6.0 :

Installation link : https://spinalcordtoolbox.com/user_section/installation.html

# Dataframe generation

To generate the dataframe run :

~~~
python generate_dataframe.py --data /path/to/CanProCo --lesion /path/to/lesion/segmentation --discs /path/to/discs/segmentation --spinal-cord /path/to/spinal/cord/segmentation --timepoint M0 --exclude-file /path/to/exclude/file --output /path/to/output/folder
~~~

> **Note**
> It uses the file `image.py` to fix the orientation of masks if they are not the same as the image's orientation

# Dataframe analysis

The following Notebook `dataframe_analysis.ipynb` details the analysis of the generated dataframe. To use it, update the link to the CSV dataframes. 
