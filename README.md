# Automated quantification of photoreceptor outer segments in developing and degenerating retinas on microscopy images across scales
This repository contains the code which was used to quantify photoreceptor outer segments (POS) from light microscopy (QuaPOS-LM) and transmission electron microscopy (QuaPOS-TEM) images as published in [link](). The read me file is separated into two parts describing each method separately.

## Light microscopy analysis (QuaPOS-LM)
Quantification of POS number, size, shape, and intensity from cryosections stained with the cone marker S-opsin and recorded with fluorescent light microscopy.

### 1 Overview
The provided random forest classifier (QuaPOS-LM) intends to predict cone POS from cryosections stained with S-opsin. Here, analysis was carried out in `Python 3.9`. The code is separated in 3 main parts. Part 1, was used to train and validate the random forest classifier. Part 2 was used to extract number, size, shape, and intensity from two datasets stained for cone POS. Part 3 contains python scripts used for statistical data analysis for the two datasets. First, postnatal development was analysed in wildtype (C57BL/6JRj, WT) mice. Second, cone POS loss was analysed in a genetic mutant mouse line, cone photoreceptor function loss 1 (Cpfl1) and compared to age-matched WT-control animals.

- authors: Florian Salomon, Robert Haase
- Institute: Center for Regenerative Therapies Dresden (CRTD), Technische Universität Dresden, Dresden, Germany

### 2 Repository contents for light microscopy analysis using QuaPOS-LM
The python scripts are stored in the folder `quapos-lm`. The folder contains the following sub-structure.

- Folder `01-training-and-validation`:
    - python scripts which were used to train and validate the random forest classifier
    - contains the random forest classifier `quapos_lm.cl`
- Folder `02-feature-extraction`:
    - subdivided into `01-wt`, and `02-cpfl`
    - the folder contains python scripts which were used to extract features from the corresponding image datasets
    - data tables are saved under folder measurements 
- Folder `03-plots-and-statistics`:
    - subdivided into `01-wt`, and `02-cpfl`
    - python scripts for the statistical analysis
    - computed plots are saved under plots-and-images respectively
- Folder `04-microscopy-images`:
    - subdivided into `01-wt`, and `02-cpfl`
    - python scripts to process images and create masks and images shown in the corresponding figures in the publication
    - respective images are saved under the folder plots-and-images respectively
 
Additionally, a python file with custom function used to process images, plots, and statistics (`quapos-lm.py`) and a yml-file (`quapos-lm.yml`) are provided.

The folder data is empty in the github repository. Respective images are provided under the following repository:  ADD LINK HERE
 
Usually, the python scripts are commented well and self-explanatory. However, if you do encounter a problem or have questions please submit a post in the issues section. 

### 3 Software training and validation QuaPOS-LM
The pixel classifier (QuaPOS-LM) was trained using [apoc](https://github.com/haesleinhuepf/apoc) (0.12.0), ground truth annotation were created using [napari](https://github.com/Napari/napari) (0.4.17). Annotation were provided for a postnatal development series of WT mice retinal sections stained for S-opsin. To train a pixel classifier which distinguishes between S-opsin signal and its background ground truth annotation were drawn. 100 pixel for each class were annotated in 3 biological replicates from 7 different developmental stages (21 images, from P08 to P24). After establishing the random forest pixel classifier, its performance was estimated. A set of ground truth images were annotated. The images were predicted with QuaPOS-LM and a confusion matrix computed by comparing the prediction and the correspodning ground truth annotation. A confusion matrix was then used to calculate different performance scores.

### 4 Workflow for light microscopy analysis
The following workflow should be applied when a new dataset should be analysed:

    1. Separate multi-channel images and save them as `.tif` file format ([fiji/imageJ software](https://imagej.net/software/fiji/downloads) is recommended since automated [imagej macro scripts](https://forum.image.sc/t/macro-in-batch-processing-to-split-channel-and-save-single-channel-image/26426/2) are available)
    2. 

### 5 System requirements
Provided notebooks and code was written and heavily tested in Windows 10 and Python 3.9. A virtual environment containing devbio-napari (0.8.1) was created using mamba (1.1.0). Some of the provided scripts rely on packages ('pyclesperanto.prototype' and 'APOC') which require a graphics card for better performance. The pixel classifier was trained and tested using APOC (0.12.0). Provided images were provided as tif file format containing a single channel of interest.

### 6 Installation guide
To use the provided code it is recommended to create a virtual environment with a conda distribution. We recommend the distribution [mamba/miniforge](https://github.com/conda-forge/miniforge#mambaforge) and using the following description of setting up mamba for your local machine as [here](https://haesleinhuepf.github.io/BioImageAnalysisNotebooks/01_introduction/readme.html). Additionally, it is recommended to install the devbio-napari environment along with seaborn `mamba create --name devbio-napari-env python=3.9 devbio-napari seaborn -c conda-forge`. Alternatively, the environment could also be recreated with the provided yml file and the command `mamba env create -f quapos-lm.yml`

### 7 Data access and how to use the code
Image data was uploaded to the following [repository](). The folder quapos-lm includes the folder `01-original-data` which contains image data in czi (Zen) file format which are unprocessed images as acquired by the microscope (Apotome ImagerZ2, 20x air objective). 
 
The classifier file `quapos-lm.cl` is stored under the folder `01-training-and-validation`. Notebooks which require the usage of the model access the classifier with a relative file path (`quapos_lm = apoc.ObjectSegmenter(opencl_filename="../../01-training-and-validation/quapos-lm.cl"`), depending on your datastructure the relative path might need to be adjusted. If retraining of the classifier is required, (e.g., if the classifier file is unable to handle data from your microscope) the folder `01-training-and-validation` contains different notebooks which make it easy to follow the original workflow and can be adapted to your dataset.

The folder `02-feature-extraction` contains jupyter notebooks which were used to extract features using [napari-simpleitk-image-processing](https://github.com/haesleinhuepf/napari-simpleitk-image-processing) (0.4.5) as well as porespy [porespy](https://github.com/PMEAL/porespy) (2.3.0). The jupyter-notebooks save the data tables as `csv` file format in the folder measurements. 

## Quantification of photoreceptor outer segment membrane stack alignment and morphology on transmission electron microscopy images (QuaPOS-TEM)

### 1 Overview
A custom MATLAB code for quantification of POS membrane stack morphology was developed on the basis of orientation and coherency analysis… (references). The code extracted the orientation of membranes and their coherency from the image gradient of sliding image patches on transmission electron microscopy (TEM) images.
TEM images of retinal sections acquired at 5000x magnification (287.5 pixels = 1 µm, acquired at FEI Morgagni D268 (camera: MegaView III, Olympus) running at 80kV acceleration voltage) were used to analyse the POS ultrastructure. Single POS were selected as separate ROIs using the Selection Brush Tool in ImageJ (version 1.54b). POS were identified as subcellular structures with electron-dense membranes at the tip of a connecting cilium or between the mitochondria-rich inner segments and the highly pigmented RPE. All ROIs of one image were analysed individually but saved together for per-specimen analysis. 

During analysis with QuaPOS-TEMthe alignment of the POS membrane stacks was calculated both locally (within a given radius) and globally (across individual POS). The coherency analysis was applied to quantify the POS morphology of WT and two inherited retinal degeneration (retinal degeneration 19 (rd19; mutation in prominin1) and rhodopsin knock-out (RhoKO)) mouse lines at 1 month of age.

Briefly, the code and attached functions contain the following steps:

Opening of files, creation of folder for processed data, reading of metadata
Definition of box radius for coherency analysis (here box radius for local coherency = 12 pixels), factors for downsampling
Extraction of ROIs into zip folder
Orientation analysis based on method gradient, saving of results file 1
Coherency analysis, saving of results file 2 (+deletion of results file 1)
Coherency analysis per ROI, saving of results file 3 (+deletion of results file 2)
Generation and saving of plots (whole image with orientation field and ROIs; whole image with coherency field, ROIs and their respective global coherency; image per ROI with orientation field; image per ROI with coherency field and global coherenc; density function of local coherency per ROI; polar histogram of coherencies)
Export of mean local coherency, global coherency and angle of global coherency into txt file
- authors: Karl Hoffmann, Suse Seidemann
- Institutes: Center for Regenerative Therapies Dresden (CRTD), Technische Universität Dresden, Dresden, Germany, MPI-CBG, Dresden, Germany

### 2 Repository Contents

- "Hauptcode"
- aufgerufene Codes/Funktionen (orientierungs und coherency analysis)
- Code for analysis of alignment of angle of global coherency per biological replicate
- used TEM images and respective files with ROIs: repository link
- retrieved files with collection of analysis data and generated plots: repository link 

### 3 System Requirements

- ImageJ (1.54b)
- MATLAB (R2023b)
- RAM (recommendation) or maybe access to computational cluster
- enough storage space for data (recommendation: external harddrive)
- Graph Pad Prism, napari or R for statistical analyses
  
### 5 Data access and workflow for using QuaPOS-TEM analysis

- install ImageJ and MATLAB
- download all MATLAB codes into one folder
- download the example data (images plus ROIs) 
  or
- use any similar TEM images (the code was optimized for a resolution of 287.5 pixels = 1 µm (5000x magnification at FEI Morgagni D268 (camera: MegaView III, Olympus) running at 80kV acceleration voltage) collected in one folder with a seperate subfolders for each biological sample
- select POS as ROIs using ImageJ and save imagename_ROIs.roi" files for each image in same folder
- adapt box radius for local alignment or dsG/dsO according to your resolution
- adapt MATLAB code to open respective folders and create output-folders
- run code
- do statistical analysis


## Citation
- How to cite the code

## Acknowledgements
- comes later
