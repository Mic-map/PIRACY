import argparse
import pathlib
from PIRACY import PIRACY

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--datapath", type=pathlib.Path, required=True)
    parser.add_argument("--subject", type=str, required=True)
    parser.add_argument("--ses", type=int, required=True)
    parser.add_argument("--run",  type=int, required=True)
    parser.add_argument("--epiPrefix",  default='task-rest_acq-EPI', type=str, required=False)
    parser.add_argument("--epiSuffix",  default='bold.nii.gz', type=str, required=False)
    parser.add_argument("--anatSuffix",  default='acq-FSEMS.nii.gz', type=str, required=False)

    args = parser.parse_args()

    print("-----------------------------------------------------------")
    myPiracy = PIRACY(args.datapath, args.subject, args.ses, args.run)
    #anat
    myPiracy.anat_scale_pixdim_x10()
    myPiracy.anat_brain_masking()
    myPiracy.anat_reg_to_atlas()
    #fMRI
    myPiracy.fmri_scale_pixdim_x10()    
    myPiracy.fmri_loose_mask()
    myPiracy.fmri_mpdenoising()
    myPiracy.fmri_topup()
    myPiracy.fmri_tight_mask()
    myPiracy.fmri_reg_to_anat()
    myPiracy.fmri_slicetiming_correction_spatial_smoothing()
    myPiracy.fmri_ICA(components=40, ndelete=10)
    myPiracy.create_connectome(gsr=True)