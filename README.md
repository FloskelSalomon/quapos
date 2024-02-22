# Automated quantification of photoreceptor outer segments in developing and degenerating retinas on microscopyic images across scales
The repository contains the code which was used to quantify outer segments in the afformentioned paper.

## 1. Overview
The provided code analyses outer segments from light microscopy and electron microscopy images. 

### Light Microscopy Analysis (QuaPOS-LM)
The provided python code intends to analyse outer segments (OS) from fluorescent light microscopy images. Analysis was carried out in Python 3.9. The provided python code generates a random forest classifier from S-opsin stained retinal sections and enables their analysis. Afterwards the random forest classifier is used to extract different features to analyse two datasets. First, postnatal development of outer segments is analysed. Second, OS degeneration is analysed in a mutant mouse model called cone photoreceptor function loss 1 (Cpfl1). The provided python codes extract different shape, size, and intensity features and computes different plots and statistics to analyse their change over time and genotype.

- authors: Florian Salomon, Robert Haase
- Institute: Center for Regenerative Therapies Dresden (CRTD), Technische Universität Dresden, Dresden, Germany

### Electron Microscopy Analysis
A custom MATLAB code for quantification of outer segment membrande stack morphology was developed on the basis of orientation and coherency analysis… (references). The code extracts the orientation of membranes and their coherency from the image gradient of sliding image patches on transmission electron microscopic images.
The alignment of the membrane stacks is calculated both locally (within a given radius) and globally (across individual OSs). The coherency analysis was applied to quantify the OS morphology of wild type and two inherited retinal degeneration (retinal degeneration 19 (rd19; mutation in prominin1) and rhodopsin knock-out (RhoKO)) mouse lines at 1 month of age.

Briefly, the code (and attached functions) contains the following steps:

Opening of files, creation of folder for processed data, reading of metadate
Definition of box radius for gradient and coherency analysis, DsO, DsG
Extraction of ROIs into zip folder
Orientation analysis based on method gradient, saving of results file 1
Coherency analysis, saving of results file 2 (+deletion of results file 1)
Coherency analysis per ROI, saving of results file 3 (+deletion of results file 2)
Generation of images (whole image with orientation field and ROIs, whole image (“gradient version”?) with coherency field, ROIS and their global coherency, image per ROI with orientation field, other type of image per ROI (gradient??) with coherency field, density functions of orientation fields per ROI)
Export of mean local coherency, global coherency and global coherency angle into txt file

- authors: Suse Seidemann, Karl Hofmann
- Institutes: Center for Regenerative Therapies Dresden (CRTD), Technische Universität Dresden, Dresden, Germany, MPI-CBG, Dresden, Germany

## 2. Repository Contents
The repository provides code for the analysis of light microscopy images with python and code for the analysis of electron microscopy images with matlab.

### Light Microscopy
The different subfolders have the following contents: 

- Folder `01-training-and-validation`:
    - python scripts used to train and validate a random forest classifier
    - contains the random forest classifier `quapos_lm.cl`
- Folder `02-feature-extraction`:
      - subdivided into 01-wt, and 02-cpfl
      - python scripts for the feature extraction and their respective post-processing
      - data tables are saved under folder measurements 
- Folder `03-plots-and-statistics` and folder `04-microscopy-images` (subdivided into 01-wt, and 02-cpfl respectively) contain different scripts which were used to quantify measured features upon extraction and create the respective plots and images. According notebooks generate plots, heatmaps, and microscopy masks which are saved respectively under the folder plots-and-images.

- `classifier.cl`: the trained supervised machine learning model to segment images into signal and background
- `model-validation.ipynb`: example workflow how the model has been validated
- `image-segmentation.ipynb`: demonstrating segmentation of the model of selected images
- `feature-extraction.ipynb`: example workflow how features were extracted from one image
- `correlation-matrix-postnatal-development.ipynb`: notebook to demonstrate how the correlation matrix for the postnatal development was computed
- `umap-postnatal.development.ipynb`:
- `statistical-analysis.ipynb`: statistical analysis of selected features shown in the paper

Additionally, the folder contains subfolders were computed images, plots, and extracted features are saved accordingly. The repository here does not contain the original czi files. To access them please use the following link: ....

### Electron Microscopy

- Beispielbilder und Beispiel-Ergebnis-files
- "Hauptcode"
- aufgerufene Codes/Funktionen (orientierungs und coherency analysis)

## 3. System Requirements

### Light Microscopy 
- python notebooks are support for Windows
- python codes were written in Windows 10
- python codes were tested for Windows 10, mamba version ..., and devbio napari version ...
- the pixel classifier was written and trained with APOC version ...
- following packages and code were used to execute the code (ask Robert whether it is possible to execute a python command, which shows the current loaded packages and their version)
    - in every notebook the current package versions used could be implemented

### Electron Microscopy
- ImageJ
- MATLAB
- RAM (recommendation) or maybe access to computational cluster
- enough storage space for data (recommendation: external harddrive)
- Graph Pad Prism, napari or R for statistical analyses

## 4. Installation Guide

### Light Microscopy
- install mamba and the devbio napari environment (add link to Roberts notebook and the yml file)

### Electron Microscopy
- install ImageJ
- install MATLAB
  
## 5. Data access + using the code
- all microscopy images were uploaded to the following repository

### Light Microscopy
- to use the python codes provided here with the light microscopy data please download the data as well as the python notebooks
- after putting them both into corresponding subfolders (one notebooks, and one data) it should be possible to run the python notebooks smoothly

### Electron Microscopy
- download all MATLAB codes into one folder
- use the example data (images plus ROIs) and dowload it into a folder
  or
- use any similar TEM images (code optimized for 287.5 pixels = 1 µm (5000x magnification at FEI Morgagni D268 (camera: MegaView III, Olympus) or a Jeol JEM1400 Plus (camera: Ruby, JEOL) both running at 80kV acceleration voltage) collected in one folder
- select ROIs using ImageJ and save imagename_ROIs.roi" files for each image in same folder
  - adapt box radius for local alignment or dsG/dsO accordingly
- adapt MATLAB code to open write folders and create output-folders
- run code
- do statistical analysis

## 6. Software Training for Pixel Classification for Analysis of Light Microscopy Images
- Training with Napari version ... and APOC version ...
- the annotation used for the software training are provided in the corresponding folder in the repository 
- annotation were drawn with the labels tool in napari version ...
- around 100 pixels in each class (signal and background) were annotated in each image
- images from different timepoints of postnatal development were used for the training
      - 7 timepoints between P8-P24
      - 3 images per timepoint

## 7. Workflow for Light Microscopy analysis
- provide the workflow with inbetween steps, can also be done in notebook

## 8. Workflow for Electron Microscopy analysis
- acquire 5000x images or similar
- select ROIs using ImageJ and save imagename_ROIs.roi" files for each image in same folder
  - adapt box radius for local alignment or dsG/dsO accordingly
- adapt MATLAB code to open write folders and create output-folders
- run code
- do statistical analysis

## 9. Citation
- How to cite the code

## 10. Acknowledgements
- comes later
