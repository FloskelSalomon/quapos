# Automated quantification of photoreceptor outer segments in developing and degenerating retinas on microscopy images across scales
This repository contains the code which was used to quantify photoreceptor outer segments from light microscopy images (QuaPOS-LM) and transmission electron microscopy images (QuaPOS-TEM) as published in [link](). This read me file is separated via the used method.

## Light microscopy analaysis (QuaPOS-LM)

### 1 Overview
The provided codes aim to analyse outer segments from light microscopy and transmission electron microscopy images.

The provided python code intends to analyse outer segments (OS) from fluorescent light microscopy images. Analysis was carried out in Python 3.9. The provided python code generates a random forest classifier from S-opsin stained retinal sections and enables their analysis. Afterwards the random forest classifier is used to extract different features to analyse two datasets. First, postnatal development of outer segments is analysed in C57BL/6JRj (wildtype) animals. Second, OS degeneration is analysed in a mutant mouse model called cone photoreceptor function loss 1 (Cpfl1) and compared to age-matched WT animals. The provided python codes extract different shape, size, and intensity features and compute different plots and statistics to analyse their change over time and genotype.

- authors: Florian Salomon, Robert Haase
- Institute: Center for Regenerative Therapies Dresden (CRTD), Technische Universität Dresden, Dresden, Germany

### 2 Repository Contents
The repository provides jupyter notebooks with python code for the analysis of S-opsin stained light microscopy images. The main folder is `quapos-lm`, it contains the different subfolders:

- Folder `01-training-and-validation`:
    - python scripts used to train and validate a random forest classifier
    - contains the random forest classifier `quapos_lm.cl`
- Folder `02-feature-extraction`:
    - subdivided into `01-wt`, and `02-cpfl`
    - python scripts for the feature extraction and their respective post-processing
    - data tables are saved under folder measurements 
- Folder `03-plots-and-statistics`:
    - subdivided into `01-wt`, and `02-cpfl`
    - python scripts for the statistical analysis
    - computed plots are saved under plots-and-images respectively
- Folder `04-microscopy-images`:
    - subdivided into `01-wt`, and `02-cpfl`
    - python scripts to postprocess image masks for publication
    - respective images saved under the folder plots-and-images respectively
 
Additionally, a python file with custom function used to process images, plots, and statistics (`quapos-lm.py`) and a file of the virtual environment (`quapos-lm.yml`) are provided.

The folder data is empty in the github repository. Respective images are provided under the following repository:  ADD LINK HERE

### 3 System requirements
Provided notebooks and code was written and heavily tested in Windows 10 and Python 3.9. A virtual environment containing devbio-napari (0.8.1) was created using mamba (1.1.0). Some of the provided scripts rely on packages ('pyclesperanto.prototype' and 'APOC') which require a graphics card for better performance. The pixel classifier was trained and tested using APOC (0.12.0).

### 4 Installation guide
To use the provided code it is recommended to create a virtual environment with a conda distribution. We recommend the distribution [mamba/miniforge](https://github.com/conda-forge/miniforge#mambaforge) and using the following description of setting up mamba for your local machine as [here](https://haesleinhuepf.github.io/BioImageAnalysisNotebooks/01_introduction/readme.html). Additionally, it is recommended to install the devbio-napari environment along with seaborn `mamba create --name devbio-napari-env python=3.9 devbio-napari seaborn -c conda-forge`. Alternatively, the environment could also be recreated with the provided yml file and the command `mamba env create -f quapos-lm.yml`

### 5 Data access and using the code
Image data was uploaded to the following [repository](). This includes the raw czi images along with the respective images and annotation (in tif format) which were used to train and validate the random forest classifier. To use the python code please download the folder `quapos-lm` in a desired location. If you would like to work with the original image dataset, download the images from its repository in the folder `data`.

### 6 Software training for the random forest pixel classifier
The pixel classifier was trained apoc (0.12.0), provided mask were created using napari (0.4.17). Annotation were provided for a postnatal development series of mice retinal sections stained for S-opsin. To train a pixel classifier which can distinguish between 2 classes (signal and background) provided images were annotated. 100 pixel for each class were annotated in 3 biological replicates from 7 different timepoints (21 images). Afterwards establishing the random forest pixel classifier, its performance was estimated. A set of ground truth images were annotated. The raw images were predicted with QuaPOS-LM and a confusion matrix computed by comparing the prediction and its annotation. The confusion matrix was then used to calculate different performance scores. 

### 7 Workflow for light microscopy analysis
A bioinformatic workflow for the light microscopy analysis can be found under folder 02-feature-extraction-workflow. It is recommended to provide files in tif format containing 1 channel. The provided code normalises the data in accordance to the training data. The normalised image will then be segmented using the pixel classifier. Finally, some features can be extracted which can be processed subsequently. The according folder 01-wt and 02-cpfl contain python scripts which require folders of image data to extract features. For the analysis here functions from `napari-simpleitk-image-processing` (0.4.5) and `porespy` (2.3.0) to extract features.

Finally, after feature extraction the folder 03-plots-and-statistics is used to analyse the data.

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
