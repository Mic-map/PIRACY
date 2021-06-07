# PIRACY
An optimized pipeline for functional connectivity analysis in the rat brain (“PIRACY”) is an image-processing pipeline for functional MRI data that has been developed to identify and correct various specific artifacts and confounding factors for an improved accuracy, precision, and robustness in functional connectivity analysis for rat. The pipeline includes:
  1. Pixdim scaling in the nifti header
  2. Brain masking
  3. MP-PCA denoising
  4. Field distortion correction (Topup)
  5. co-registration
  6. Brain parcellation
  7. Slice-timing correction and spatial smoothing
  8. Melodic ICA
  9. FIX auto-classification and cleaning
  10. Functional connectivity creation

PIRACY is modular and can easily be tuned to study-specific needs or requirements. The current implementation is written in Python, but requires the installation of Matlab, SPM, ANTs and FSL. The following MATLAB toolboxes are also required to install by users: SPM, MP-PCA denoising (https://github.com/NYU-DiffusionMRI/mppca_denoise) and NIFTI toolbox (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image). To call MATLAB in python, Matlab's python engine should be installed (https://www.mathworks.com/help/matlab/matlab-engine-for-python.html).

In order to enable FIX ICA auto-classification, specific configuration is required. Please refer to FIX's userguide for details (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX/UserGuide). 

The PIRACY pipeline has been evaluated and described in more technical detail in following publication:
Yujian Diao, Ting Yin, Rolf Gruetter and Ileana O Jelescu, PIRACY: An Optimized Pipeline for Functional Connectivity Analysis in the Rat Brain. 
Front. Neurosci. 15, (2021). (https://doi.org/10.3389/fnins.2021.602170)

# Recommended Usage
It is recommended to run PIRACY in a python script following the example in "test.py" although it is possible to run in a command line, e.g,
    python3 /yourpath/PIRACY/main.py --step='preprocessing' --datapath=/your/data/path --subject='CTL3693' --ses=1 --run=1
    
We assume the data is named and stored following the rules of Brain Imaging Data Structure (BIDS, https://bids.neuroimaging.io/) 
For example, for the data of subject CTL3692:

![image](https://user-images.githubusercontent.com/38806138/114040456-ed0b2380-9883-11eb-9700-bc804ba1fe95.png)

