# Automated quantification of photoreceptor outer segments in developing and degenerating retinas on microscopy images across scales
This repository contains the code which was used to quantify photoreceptor outer segments (POS) from light microscopy (QuaPOS-LM) and transmission electron microscopy (QuaPOS-TEM). The corresponding image datasets are distirbuted under the BioImage Archive and can be accessed [here](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1078). 

Here, the two methods are described separately from each other.

## Quantification of photoreceptor outer segment (POS) number, size, shape, and intensity from light microscopy images (QuaPOS-LM)
Quantification of POS number, size, shape, and intensity from cryosections stained with the cone marker S-opsin and recorded with fluorescent light microscopy.

### 1 Overview
The provided random forest classifier (QuaPOS-LM) intends to predict cone POS from cryosections stained with S-opsin. Here, analysis was carried out in `Python 3.9.16`. The code is separated in 3 main parts. Part 1, was used to train and validate the random forest classifier. Part 2 was used to extract number, size, shape, and intensity from two datasets stained for cone POS. Part 3 contains python scripts used for statistical data analysis for the two datasets. First, postnatal development was analysed in wildtype (C57BL/6JRj, WT) mice. Second, cone POS loss was analysed in a genetic mutant mouse line, cone photoreceptor function loss 1 (Cpfl1) and compared to age-matched WT-control animals.

- authors: Florian Salomon $^1$, Robert Haase $^2$
- Institute:
    - 1: Center for Regenerative Therapies Dresden (CRTD), Technische Universität Dresden, Dresden, Germany
    - 2: Center for scalable data analytics and artificial intelligence, Universität Leipzig, Leipzig, Germany

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
The pixel classifier (QuaPOS-LM) was trained using [apoc](https://github.com/haesleinhuepf/apoc) (0.12.0), ground truth annotation were created using [napari](https://github.com/Napari/napari) (0.4.17). The classifier was trained with the following settings:
- number of trees: 100
- number of decisions: 2
- filter operations:
    - gaussian blur (sigma 1)
    - difference of gaussian blur (sigma 1)
    - laplace box of gaussian blur (sigma 1)

Sparse annotation were provided for a set of images from a postnatal development series of WT mice retinal sections stained for S-opsin (3 biological replicates from 7 developmental stages, 21 images). 100 pixel for each class (signal and background) were annotated in each image. After establishing the binary classifier, its performance was estimated. A set of ground truth images were annotated. The images were predicted with QuaPOS-LM and a confusion matrix computed by comparing the prediction and the correspodning ground truth annotation. A confusion matrix was then used to calculate different performance scores.

### 4 Workflow for light microscopy analysis
The following workflow was used to analyse new image datasets and is recommended when a new dataset should be analysed:

- Separate multi-channel images and save them as `.tif` file format ([fiji/imageJ software](https://imagej.net/software/fiji/downloads) is recommended since automated [imagej macro scripts](https://forum.image.sc/t/macro-in-batch-processing-to-split-channel-and-save-single-channel-image/26426/2) are available).
    - Note: If several technical replicates per biological replicate are acquired it is recommended to add this information to the file name (e.g., all images of animal 13 get the tag `-biological-replicate-13` in their filename). This will make it possible to add a column with the corresponding biological replicate to the dataframe and group the data by that column later.
- Store image data for analysis in a folder
- Extract feauters, to extract features from unprocessed images QuaPOS-LM provides the following workflow:
    - background subtraction (top-hat-filter)
    - intensity normalisation
    - label prediction
    - extract features
- Process the dataset
    - Rescale size features
    - Filter the dataset (remove inf values and filter 5- and 95-percentiles)
    - Calculate average values (if applicable, for biological and technical replicates)
    - Add additional features (e.g., number of POS labels, and summed POS volume)
- Statistical data analysis

### 5 System requirements
The provided code was written in Python (3.9.16) [jupyter-lab](https://github.com/jupyterlab/jupyterlab) (3.6.1). The provided notebooks were heavily tested in Windows 10. The pixel classifier was trained and applied to image data using [APOC](https://github.com/haesleinhuepf/apoc). The notebooks rely on different software packagaes from [devbio-napari](https://github.com/haesleinhuepf/devbio-napari) (version 0.8.1 was used here) which was installed in a virtual environment using [mamba/miniforge](https://github.com/conda-forge/miniforge#mambaforge) (1.1.0). Additionally, [seaborn](https://github.com/mwaskom/seaborn) (> 0.13.2) is required to compute various plots. A RAM of at least 8 GB is recommended. Some packages could require a graphics card, if you run into any troubles please refer to [devbio-napari](https://github.com/haesleinhuepf/devbio-napari).

### 6 Installation guide
You can install all softwar packages using conda/mamba. If you have not used conda before, please read [this guide first](https://biapol.github.io/blog/mara_lampert/getting_started_with_mambaforge_and_python/readme.html). Afterwards, you can create a new environment by entering the following command into your prompt:

    mamba create --name pos-analysis-quapos python=3.9 devbio-napari seaborn -c conda-forge

If the installation was successfull you can load the environment with the following command:

    mamba activate pos-analaysis-quapos

And start jupyter-lab with the following command:

    jupyter-lab

Alternatively, the environment could also be recreated using the provided `.yml` file in the folder `quapos-lm`.

### 7 Data access and how to use the code
To use the provided code please download the folder `quapos-lm` and store it on your local machine. Afterwards you should be able to find it under the stored location using the file browser in the jupyter lab interface. The light microscopy image dataset analysed here was uploaded to the following [repository](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1078) under quapos-lm. Please download the data and provide it in the folder `data`. The folder quapos-lm includes the folder `01-original-data` which contains image data in czi (Zen) file format which are unprocessed images as acquired by the microscope (Apotome ImagerZ2, 20x air objective) as well as the folder `02-data-for-pixel-classifier` where single channels and corresponding datasets used for testing and training the classifier are stored.
 
The classifier file `quapos-lm.cl` is stored under the folder `01-training-and-validation`. Notebooks which require the usage of the model access the classifier with a relative file path (`quapos_lm = apoc.ObjectSegmenter(opencl_filename="../../01-training-and-validation/quapos-lm.cl"`), depending on your datastructure the relative path might need to be adjusted. If retraining of the classifier is required, (e.g., if the classifier file is unable to handle data from your microscope) the folder `01-training-and-validation` contains different notebooks which make it easy to follow the original workflow and can be adapted to your dataset.

The folder `02-feature-extraction` contains jupyter notebooks which were used to extract features using [napari-simpleitk-image-processing](https://github.com/haesleinhuepf/napari-simpleitk-image-processing) (0.4.5) as well as porespy [porespy](https://github.com/PMEAL/porespy) (2.3.0). The jupyter-notebooks save the data tables as `csv` file format in the folder measurements. Additionally, the folder contains also different jupyter notebooks which were used to process the dataset. To provide your own images to the classifier, it should only be required to change the directories and filenames of the dataset in the according line of the jupyter notebook.

After the feature extraction and preprocessing of the data tables of your own data you can start with the statistical data analysis (`03-plots-and-statistics`). The provided notebooks show plots created with [matplotlib](https://github.com/matplotlib/matplotlib) (3.7.0), and [seaborn](https://github.com/mwaskom/seaborn) (> 0.13.2) as well as corresponding statistical analysis using [scipy](https://github.com/scipy/scipy) (1.10.1) and [statsmodels](https://github.com/statsmodels/statsmodels) (0.14.1). In the analysis shown in the paper, only selected features were analysed. If you are interested in investigating something else you can follow the notebooks in the provided folders. Additionally, the folder `02-cpfl` provides notebooks for the analysis of two different groups over time, and can be adjusted accordingly to the groups of your own experiments. 

## Quantification of photoreceptor outer segment membrane stack alignment and morphology on transmission electron microscopy images (QuaPOS-TEM)

### 1 Overview
A custom MATLAB code for quantification of POS membrane stack morphology was developed on the basis of orientation and coherency analysis. 
\[compare Karen Soans et. al. (2022) [Matrix topology guides collective cell migration in vivo](https://www.cell.com/current-biology/pdf/S0960-9822(22)01503-2.pdf). Current Biology, and Karl B. Hoffmann (2022) [Robust Identification of Topological Defects in Discrete Vector Fields with Applications to Biological Image Data (Doctoral dissertation)](https://nbn-resolving.org/urn:nbn:de:bsz:14-qucosa2-857202). Technical University Dresden\]

The code extracts the orientation of membranes in transmission electron microscopy (TEM) images as the orientation of high-intensity image features, using the image gradient in 5-by5 image patches according to the Scharr operator optimized for orientational bias \[Scharr, Hanno (2000). [Optimale Operatoren in der digitalen Bildverarbeitung (Doctoral dissertation)](https://archiv.ub.uni-heidelberg.de/volltextserver/962/1/Diss.pdf).  Universität Heidelberg\].
These orientations are nematic vectors, and their coherency (or degree of mutual alignment) is measured by the scalar nematic order parameter.

TEM images of retinal sections were acquired at 5000x magnification (287.5 pixels = 1 µm, FEI Morgagni D268 (camera: MegaView III, Olympus) running at 80kV acceleration voltage) and used to analyse the POS ultrastructure. Single POS were selected as separate ROIs using the Selection Brush Tool in ImageJ (version 1.54b), typicially multiple within each image. POS were identified as subcellular structures with electron-dense membranes at the tip of a connecting cilium or between the mitochondria-rich inner segments and the highly pigmented RPE. All ROIs of one image were analysed individually but saved together for per-specimen analysis. 

During analysis with QuaPOS-TEM the alignment of the POS membrane stacks was calculated both locally (within a given radius) and globally (across individual POS). The coherency analysis was applied to quantify the POS morphology of WT and two inherited retinal degeneration (retinal degeneration 19 (rd19; mutation in prominin1) and rhodopsin knock-out (RhoKO)) mouse lines at 1 month of age.

Briefly, the code and attached functions contain the following steps:

Within folders that can group images of the same condition and same age
- Process each image in the folder:
    - Definition of box radius for coherency analysis (here box radius for local coherency = 12 pixels), factors for downsampling (here: 1)
    - Orientation analysis based on method gradient across the whole image, saving of results file 1
    - Coherency analysis across the whole image, saving of results file 2 (+ deletion of results file 1)
    - Coherency analysis per ROI, calculating especially the mean local coherency of each ROI and the global coherency of each ROI, saving of results file 3 (+ deletion of results file 2)
- Generation and saving of plots
    - whole image with orientation field and ROIs
    - whole image with coherency field, ROIs and their respective global coherency
    - one image per ROI with orientation field
    - one image per ROI with coherency field and global coherenc
    - density function plots of local coherency per ROI
    - polar histogram of coherencies per ROI
- Export of mean local coherency, global coherency and angle of global coherency for all ROIs in all images into a joint file for subsequent statistical analysis
 
- authors: Karl Hoffmann $^3$, Suse Seidemann $^4$
- Institute:
    - 3: Center for Systems Biology Dresden (CSBD), Max Planck Institute of Molecular Cell Biology and Genetics (MPI-CBG), and Faculty of Computer Science, Technische Universität Dresden, all Dresden, Germany
    - 4: Center for Regenerative Therapies Dresden (CRTD), Technische Universität Dresden

### 2 Repository Contents

- main file QuaPOS_TEM.m
- functional dependencies for orientation and coherency analysis
- used TEM images and respective files with ROIs: repository link
- retrieved files with collection of analysis data and generated plots: repository link 

### 3 System Requirements

- ImageJ (1.54b)
- MATLAB (R2023b)
- sufficient RAM (depending on image size) or optionally access to computational cluster
- sufficient storage space for intermediate data (recommendation: external harddrive)
- Graph Pad Prism or other software for statistical analyses
  
### 4 Data access and workflow for using QuaPOS-TEM analysis

- install ImageJ and MATLAB
- download all MATLAB codes into one folder
- download the example data (images plus ROIs) 
  or
- use any similar TEM images (the code was optimized for a resolution of 287.5 pixels = 1 µm (5000x magnification at FEI Morgagni D268 (camera: MegaView III, Olympus) running at 80kV acceleration voltage) collected in one folder with a seperate subfolder for each biological sample
- select POS as ROIs using ImageJ and save imagename_ROIs.roi" files, keeping images of equal condition and age in one folder
- adapt MATLAB code in QuaPOS_TEM to
    -  open your respective folders and write resutls to your output folder
    -  use box radius for local alignment or dsG/dsO according to your resolution
- run QuaPOS-TEM. Intermediate result files allow for check-pointing and are re-used when the analysis is interrupted and resumed later.
- perform statistical analysis

## Issues

If you encounter any issues with the code or accessing the microscopy dataset, please file an issue in the github repository here.

## Access to Image Dataset

The image datasets which were analysed here were distributed under the BioImage Archive:

Florian Salomon, Suse Seidemann, Karl B. Hoffmann, Thomas Kurth, Ivo F. Sbalzarini, Robert Haase, Marius Ader (2024). Automated quantification of photoreceptor outer segments in developing and degenerating retinas on microscopy images across scales. BioStudies, S-BIAD1078. [https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1078](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1078) 

## Citation

If you use the code in your work, please cite the github repository as follows:

Salomon, F., Seidemann, S., Hoffmann, K. B., Kurth, T., Sbalzarini, I. F., Haase, R., & Ader, M. (2024). FloskelSalomon/quapos: Automated quantification of photoreceptor outer segments in developing and degenerating retinas on microscopy images across scales (v0.0.1). Zenodo. [https://doi.org/10.5281/zenodo.10794252](https://doi.org/10.5281/zenodo.10794252)

## License
Distributed under the terms of [BSD-3](https://opensource.org/license/BSD-3-Clause) license, "quapos" is free and open-source software.
