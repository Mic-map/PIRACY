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

    #Setup the parameters for 1 subject
    args.datapath = '/piracy_test'
    args.subject = 'CTL3693'
    args.ses = 1
    args.run = 1

    print("-----------------------------------------------------------")
    myPiracy = PIRACY(args.datapath, args.subject, args.ses, args.run)
    #over timepoints (ses)
    for ses in [1,2,3]:
        myPiracy.set_ses(ses)
        myPiracy.anat_scale_pixdim_x10()
        myPiracy.anat_brain_masking()
        #myPiracy.get_anat_brain_and_mask()
        myPiracy.anat_reg_to_atlas()
		
	#over runs
        for run in [1,2]:
            myPiracy.set_run(run)
            myPiracy.fmri_scale_pixdim_x10()
            myPiracy.fmri_loose_mask()
            myPiracy.fmri_mpdenoising()
            myPiracy.fmri_topup()
            myPiracy.fmri_tight_mask()
            myPiracy.fmri_reg_to_anat()
            myPiracy.fmri_slicetiming_correction_spatial_smoothing()
            myPiracy.fmri_ICA(components=40, ndelete=10)
	    myPiracy.fix_classify(thresholds=[20, 70])
			
	    ### FIX cleaning. 
	    ### Please manually create the noise file "fix_noise_file.txt" based on the auto-classifications of the 2 thresholds before proceeding.
	    #myPiracy.fmri_procs_done(['px', 'dn', 'topup', 'tm', 'spm']) #specify steps have been done
	    #myPiracy.fix_clean("fix_noise_file.txt")
            #myPiracy.create_connectome()    
