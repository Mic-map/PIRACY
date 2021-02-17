import argparse
import pathlib
from PIRACY import PIRACY

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--datapath", type=pathlib.Path, required=False)
    parser.add_argument("--subject", type=str, required=False)
    parser.add_argument("--ses", type=int, required=False)
    parser.add_argument("--run",  type=int, required=False)
    parser.add_argument("--epiPrefix",  default='task-rest_acq-EPI', type=str, required=False)
    parser.add_argument("--epiSuffix",  default='bold.nii.gz', type=str, required=False)
    parser.add_argument("--anatSuffix",  default='acq-FSEMS.nii.gz', type=str, required=False)

    args = parser.parse_args()

    args.datapath = '/piracy_test'
    args.subject = 'CTL3693'
    args.ses = 1
    args.run = 1

    print("-----------------------------------------------------------")
    myPiracy = PIRACY(args.datapath, args.subject, args.ses, args.run)
    for ses in [1,2,3]:
        myPiracy.set_ses(ses)
        #myPiracy.anat_scale_pixdim_x10()
        #myPiracy.anat_brain_masking()
        #myPiracy.get_anat_brain_and_mask()
        #myPiracy.anat_reg_to_atlas()

        for run in [1,2]:
            myPiracy.set_run(run)
            myPiracy.fmri_scale_pixdim_x10()
            
            # myPiracy.fmri_loose_mask()
            # myPiracy.fmri_mpdenoising()
            # myPiracy.fmri_topup()
            # #myPiracy.fmri_procs_done(['px', 'dn', 'topup', 'tm', 'spm', 'clean'])
            # myPiracy.fmri_tight_mask()
            # myPiracy.fmri_reg_to_anat()
            # myPiracy.fmri_slicetiming_correction_spatial_smoothing()
            # myPiracy.fmri_ICA(components=50, ndelete=5)
            # #myPiracy.create_connectome()