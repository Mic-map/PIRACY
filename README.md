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

PIRACY is modular and can easily be tuned to study-specific needs or requirements. The current implementation is written in Python, but requires the installation of Matlab, SPM, ANTs and FSL.

In order to enable FIX ICA auto-classification, specific configuration is required. Please refer to FIX's userguide for details (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX/UserGuide).

Matlab's python engine is also required (https://www.mathworks.com/help/matlab/matlab-engine-for-python.html).

# Recommended Usage

