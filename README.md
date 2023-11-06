# Quantification-of-Outer-Segments-on-LM-and-EM-images (will become the title of the paper later)
- repository contains the code used to analyse outer segments of photoreceptors from LM and EM images

## 1. Overview
The provided code analysis outer segments from light microscopy and electron microscopy images. 

### Light Microscopy Analysis
The provided python code intends to analyse outer segments (OS) from fluorescent ligth microscopy images. The code was applied to retinal sections from mice immunostained for the cone OS marker S-opsin and imaged with a Zeiss Apotomoe Imager Z2 and a 20x air objective. The python code makes use of a trained supervised machine learning model to segment those stained 3D images into background and signal. Afterwards, different features are measured with napari-simpleitk from segmented images and shape, intensity, and size parameters analysed.

- authors: Florian Salomon, Robert Haase
- Institute: Center for Regenerative Therapies Dresden (CRTD), Technische Universität Dresden, Dresden, Germany

According notebooks of the python code generate the following data tables:

- "features-postnatal-development.csv": features which were measured from a postnatal development series from black 6 mice between P08 and P24
- "features-wt-cpfl-comparison.csv": features which were measured from different timepoints of black 6 and cpfl1 mice

The python code generates the following plots as images:

- add later

The python code generates the following images:

- add later

### Electron Microscopy Analysis
- description by Suse
- authors: Suse Seidemann, Karl Hofmann
- Institutes: Center for Regenerative Therapies Dresden (CRTD), Technische Universität Dresden, Dresden, Germany, MPI-CBG, Dresden, Germany

## 2. Repository Contents
The repository provides code for the analysis of light microscopy images with python and code for the analysis of electron microscopy images with matlab.

### Light Microscopy
This folder contains the according notebooks to analyse S-opsin stained outer segments. The demo folder includes:

- `classifier.cl`: the trained supervised machine learning model to segment images into signal and background
- `model-validation.ipynb`: example workflow how the model has been validated
- `image-segmentation.ipynb`: demonstrating segmentation of the model of selected images
- `feature-extraction.ipynb`: example workflow how features were extracted from one image
- `correlation-matrix-postnatal-development.ipynb`: notebook to demonstrate how the correlation matrix for the postnatal development was computed
- `umap-postnatal.development.ipynb`:
- `statistical-analysis.ipynb`: statistical analysis of selected features shown in the paper

Additionally, the folder contains subfolders were computed images, plots, and extracted features are saved accordingly. The repository here does not contain the original czi files. To access them please use the following link: ....

## 3. System Requirements

### Light Microscopy 
- python notebooks are support for Windows
- python codes were written in Windows 10
- python codes were tested for Windows 10, mamba version ..., and devbio napari version ...
- the pixel classifier was written and trained with APOC version ...
- following packages and code were used to execute the code (ask Robert whether it is possible to execute a python command, which shows the current loaded packages and their version)
    - in every notebook the current package versions used could be implemented

### Electron Microscopy
- added by Suse, hopefully

## 4. Installation Guide

### Light Microscopy
- install mamba and the devbio napari environment (add link to Roberts notebook and the yml file)

### Electron Microscopy
- added by Suse later

## 5. Data access + using the code
- all microscopy images were uploaded to the following repository

### Light Microscopy
- to use the python codes provided here with the light microscopy data please download the data as well as the python notebooks
- after putting them both into corresponding subfolders (one notebooks, and one data) it should be possible to run the python notebooks smoothly

### Electron Microscopy
- added by Suse, again

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
- added by Suse

## 9. Citation
- How to cite the code
