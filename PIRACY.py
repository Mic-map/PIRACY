
import matlab.engine
import os
import glob
import nibabel as nib
from tqdm import tqdm
from collections import OrderedDict
from functools import reduce
from utils import (execute_cmd, fmri_masking, adjust_mask, fmri_topup, ants_reg_fmri_to_anat,
                         ants_reg_anat_to_atlas, fmri_melodic_ica, create_dummy_mcf_par,fix_ica_classify, fmri_ica_clean_noise)


eng = matlab.engine.start_matlab() 

class PIRACY(object):
    def __init__(self, datapath, sub, ses, run, fmri_prefix = 'task-rest_acq-EPI', fmri_suffix='bold.nii.gz', 
                    fmri_revPE_marker='revPE', anat_suffix='acq-FSEMS.nii.gz', f_anat='anat', f_fmri = 'func'):
        self.datapath = datapath
        self.f_anat = f_anat
        self.f_fmri = f_fmri
        self.subject = sub
        self.subject_type = sub[:3]         
        self.ses = ses
        self.run = run 
        self.fmri_prefix = fmri_prefix
        self.fmri_suffix = fmri_suffix
        self.anat_suffix = anat_suffix
        self.fmri_revPE_marker = fmri_revPE_marker
        self.anat_brain_mask = ''
        self.anat_brain = ''
        self.fmri_brain_loosemask = ''
        self.fmri_brain_tightmask = ''
        self.fmri_proc_suffix_dict = OrderedDict([
            ('px', '_x10'),    #pixdim scaling by 10
            ('dn', '_dn'),     #denoising
            ('topup', '_uw'),  #topup
            ('tm', '_m'),      #tight masking
            ('spm', '_sa'),    #spm: slice timing correction and spatial smoothing
            ('clean', '_cl'),  #ICA cleaning
            ('reg', '_label')  #fMRI registering to anat, mapping back ROI labels

        ])

        self.set_file_names_and_dirs()
        self.anat_roi_labels = ''
        self.fmri_roi_labels = ''

    def set_file_names_and_dirs(self):
        #self.proc_suffix=''
        self.f_sub = "sub-{}".format(self.subject)
        self.f_ses = "ses-{}".format(self.ses)
        self.run_ = "run-{}".format(self.run)
        self.fmri_name = "{0}_{1}_{2}_{3}_{4}".format(self.f_sub, self.f_ses, self.fmri_prefix, self.run_, self.fmri_suffix)
        self.fmri_revPE_name = "{0}_{1}_{2}_{3}_{4}_{5}".format(self.f_sub, self.f_ses, self.fmri_prefix, self.run_, self.fmri_revPE_marker, self.fmri_suffix)
        self.anat_name = "{0}_{1}_{2}".format(self.f_sub, self.f_ses, self.anat_suffix)
        self.fmri_dir = reduce(os.path.join, [self.datapath, self.f_sub, self.f_ses, self.f_fmri])
        self.anat_dir = reduce(os.path.join, [self.datapath, self.f_sub, self.f_ses, self.f_anat])
        self.fmri_proc_name = self.fmri_name
        self.fmri_proc_revPE_name = self.fmri_revPE_name

    def fmri_scale_pixdim_x10(self, suffix = '_x10'):
        print('--------1.Scale the pixel size in header of fmri data (px)-------------')
        suffix_ = "{}.nii".format(suffix)
        self.fmri_proc_name  = eng.piracy_nii_scale_pixdim_x10(self.fmri_dir, self.fmri_name, suffix_)
        self.fmri_proc_revPE_name = eng.piracy_nii_scale_pixdim_x10(self.fmri_dir, self.fmri_revPE_name, suffix_)
        #self.proc_suffix = suffix

    def anat_scale_pixdim_x10(self, suffix = ''):   
        if not suffix: suffix = self.fmri_proc_suffix_dict['px']
        print('----------1.Scaling the pixel size in header of anatomical data----------')
        suffix_ = "{}.nii".format(suffix) 
        eng.piracy_nii_scale_pixdim_x10(self.anat_dir, self.anat_name, suffix_)


    def fmri_loose_mask(self, suffix='', maskSuffix='_loosemask', threshold=0.4, kernel=(5,5), iterations=3):
        print("---------2.Ceating loose mask on fMRI image ----------")
        dir_backup = os.getcwd()
        os.chdir(self.fmri_dir)
        self.fmri_brain_loosemask, _ = fmri_masking(self.fmri_proc_name, suffix=suffix, maskSuffix=maskSuffix, applyMask=False, threshold=threshold)
        #dilate mask
        adjust_mask(self.fmri_brain_loosemask, kernel, iterations)
        os.chdir(dir_backup)
        print("---Loose mask created: ", self.fmri_brain_loosemask) 

    def fmri_mpdenoising(self, suffix=''):
        print("---------3.MP-PCA denoising on fMRI images (dn)----------")
        if not suffix: suffix = self.fmri_proc_suffix_dict['dn']
        self.fmri_proc_name, self.fmri_proc_revPE_name = eng.piracy_mpdenoise(self.fmri_dir, self.fmri_proc_name, 
                                                                self.fmri_proc_revPE_name, self.fmri_brain_loosemask, suffix, True)
        #suffix_ = "{}.nii".format(suffix)
        #self.fmri_proc_name = self.fmri_proc_name.replace('.nii', suffix_)
        #self.fmri_proc_revPE_name = self.fmri_proc_revPE_name.replace('.nii', suffix_)
    def fmri_topup(self, suffix=''):
        print("---------4.Correcting field distortion of fMRI images (topup)----------")
        if not suffix: suffix = self.fmri_proc_suffix_dict['topup']
        self.fmri_proc_name = fmri_topup(self.fmri_dir, self.fmri_proc_name, self.fmri_proc_revPE_name, 
                                            suffix=suffix, acqp='acq_pars.txt', mycnf='mycnf_x10.cnf')

    def fmri_tight_mask(self, suffix='', maskSuffix='_tightmask', threshold=0.5):
        print("---------5.Creating tight brain mask on fMRI image (tm)----------")
        if not suffix: suffix = self.fmri_proc_suffix_dict['tm']
        dir_backup = os.getcwd()
        os.chdir(self.fmri_dir)
        self.fmri_brain_tightmask, self.fmri_proc_name = fmri_masking(self.fmri_proc_name, suffix=suffix, maskSuffix=maskSuffix, applyMask=True, threshold=threshold)
        os.chdir(dir_backup)    
        print("---Tight mask created: ", self.fmri_brain_tightmask)    

    def anat_brain_masking(self, prevSuffix='_x10.nii', bet_threshold=0.5, bfc=True):
        print(f"---------6.Creating brain mask on anatomical image '{self.anat_name}'------------")
        
        dir_backup = os.getcwd()
        os.chdir(self.anat_dir)
        assert os.path.exists(self.anat_name), "File not found!"
        nii_ = self.anat_name.replace(".nii", prevSuffix)
        nii_basename = nii_.split('.nii')[0]
           
        # Bias field correction
        if bfc:
            execute_cmd('fast', '-t 2 -n 3 -H 0.1 -I 4 -l 20.0 --nopve -B -o {0}'.format(nii_))
            execute_cmd('rm', '{0}'.format(nii_.replace('.nii', '_seg.nii')))
            nii_basename = "{}_restore".format(nii_basename)
        # creating mask   
        execute_cmd('bet', '{0} {0} -f {1} -g 0 -n -m'.format(nii_basename, bet_threshold))       
        # masking anatomical image
        execute_cmd("fslmaths", "{0}" \
                    " -mas {0}_mask {0}_brain".format(nii_basename)) 
        os.chdir(dir_backup)   
        self.anat_brain_mask = "{0}_mask.nii.gz".format(nii_basename)
        self.anat_brain = "{0}_brain.nii.gz".format(nii_basename) 

    def anat_reg_to_atlas(self, suffix=''):
        print("---------7.Registering anat brain to atlas template ------------")    
        if not suffix: suffix = self.fmri_proc_suffix_dict['reg']
        self.get_anat_brain_and_mask()
        self.anat_roi_labels = ants_reg_anat_to_atlas(self.anat_dir, self.anat_brain, regFolder='reg', label_suffix=suffix)

    def fmri_reg_to_anat(self, suffix='', maskSuffix='_tightmask'):
        print("---------8.Registering fmri to anat brain ------------")    
        if not suffix: suffix = self.fmri_proc_suffix_dict['reg']
        if not os.path.isfile(self.anat_brain): 
            self.get_anat_brain_and_mask()
            anat_brain_file = os.path.join(self.anat_dir, self.anat_brain)

        if not self.anat_roi_labels: 
            anat_label = anat_brain_file.replace('.nii', '_label.nii')
            assert os.path.isfile(anat_label), "----The anat ROI label file was not found: {} !!".format(anat_label)
        else:
            anat_label = os.path.join(self.anat_dir, self.anat_roi_labels)

        self.fmri_roi_labels = ants_reg_fmri_to_anat(self.fmri_dir, anat_brain_file, 
                                                        self.fmri_proc_name, anat_label, label_suffix=suffix, applyOnly=False) 
        os.chdir(self.fmri_dir)
        mask_ = glob.glob(self.fmri_name.replace('.nii', '*{}.nii'.format(maskSuffix)))
        assert mask_, "----fmri mask not found!!"
        mask = mask_[0]
        execute_cmd("fslmaths", "{0}" \
                    " -mas {1} {0}".format(self.fmri_roi_labels, mask)) 

    def fmri_slicetiming_correction_spatial_smoothing(self, suffix=''):
        print("---------9.Slice timing correction and spatial smoothing (spm)------------")
        if not suffix: suffix = self.fmri_proc_suffix_dict['spm']
        self.fmri_proc_name = eng.piracy_spm(self.fmri_dir, self.fmri_proc_name, suffix)

    def fmri_ICA(self, fix=True, components=40, ndelete=10, cutoff=100, dummyMCF=True):
        print("---------10.melodic ICA------------")
        fmri_melodic_ica(self.fmri_dir, self.fmri_proc_name, outComponents=components, 
                            useFix=fix, cutoff=cutoff, ndelete=ndelete)
        if fix:
            fmri_bsname = self.fmri_proc_name.split('.nii')[0]
            fmri_ica = os.path.join(self.fmri_dir, "{}.ica".format(fmri_bsname))
            if dummyMCF: create_dummy_mcf_par(fmri_ica)

    def fix_classify(self, thresholds=[20, 30, 70], model_file=''):
        print("---------11.1.FIX ICA noise classification------------")
        if not isinstance(thresholds, list): thresholds = [thresholds]
        fmri_bsname = self.fmri_proc_name.split('.nii')[0]
        ica = "{}.ica".format(fmri_bsname)

        for threshold in thresholds:
            print("FIX threshold: ", threshold)
            fix_ica_classify(self.fmri_dir, ica, threshold, model=model_file)

    def fix_clean(self, fix_noise_file, suffix=''):
        print("---------11.2.FIX ICA cleaning (clean)------------")
        if not suffix: suffix = self.fmri_proc_suffix_dict['clean']
        fmri_bsname = self.fmri_proc_name.split('.nii')[0]
        ica = "{}.ica".format(fmri_bsname)    
        
        self.fmri_proc_name = fmri_ica_clean_noise(self.fmri_dir, ica, noise_file=fix_noise_file, suffix=suffix)

    def create_connectome(self, label_suffix='', output_folder='connectome', gsr=True):
        print("---------12.Generate functional connectivity------------")
        os.chdir(self.fmri_dir)
        if self.fmri_roi_labels:
            fmri_label = self.fmri_roi_labels
        else:
            suffix = label_suffix if label_suffix else self.fmri_proc_suffix_dict['reg']
            fmri_label_ = glob.glob(self.fmri_name.replace('.nii', '*{}.nii'.format(suffix)))
            fmri_label = fmri_label_[0]
        eng.piracy_make_connectome(self.fmri_dir, self.fmri_proc_name, fmri_label, output_folder, gsr)

    def set_ses(self, ses):
        self.ses = ses
        self.set_file_names_and_dirs()

    def set_run(self, run):
        self.run = run
        self.set_file_names_and_dirs()

    def set_subject_type(self, subType):
        self.subject_type = subType

    def set_prefix_suffix(self, fmri_prefix='', fmri_suffix='', anat_suffix='', fmri_revPE_marker=''):
        if fmri_prefix: self.fmri_prefix = fmri_prefix
        if fmri_suffix: self.fmri_suffix = fmri_suffix
        if anat_suffix: self.anat_suffix = anat_suffix
        if fmri_revPE_marker: self.fmri_revPE_marker = fmri_revPE_marker
        if fmri_prefix or fmri_suffix or anat_suffix or fmri_revPE_marker:
            self.set_file_names_and_dirs()

    def get_anat_brain_and_mask(self, brainSuffix='_brain', maskSuffix='_mask', prevSuffix='_x10'):
        dir_backup = os.getcwd()
        os.chdir(self.anat_dir)
        anat_scaled = self.anat_name.replace('.nii', '{}.nii'.format(prevSuffix))
        assert os.path.isfile(anat_scaled), "---'{}' not found!!".format(anat_scaled)
        brain_ = glob.glob(anat_scaled.replace('.nii', '*{}.nii'.format(brainSuffix)))
        mask_ = glob.glob(anat_scaled.replace('.nii', '*{}.nii'.format(maskSuffix)))
        assert len(brain_)==1, "----Wrong/no anatomical brain image was found in '{}' !!".format(self.anat_dir)
        assert len(mask_)==1, "----Wrong/no anatomical mask was found in '{}' !!".format(self.anat_dir)
        self.anat_brain = brain_[0]
        self.anat_brain_mask = mask_[0]
        os.chdir(dir_backup) 

    def fmri_procs_done(self, procs_done=['px']):
        assert isinstance(procs_done, list)
        dir_backup = os.getcwd()
        os.chdir(self.fmri_dir)

        suffix_ = []
        for pr in procs_done:
            suffix_.append(self.fmri_proc_suffix_dict[pr])
        suffix = ''.join(suffix_)
        fmri_proc = self.fmri_name.replace('.nii', "{}.nii".format(suffix))
        if os.path.isfile(fmri_proc):
            self.fmri_proc_name = fmri_proc
        else:
            assert os.path.isfile(fmri_proc), ("---'{0}' is nonexistent in" 
                                                "{1} !!".format(fmri_proc, self.fmri_dir))
        suffix_ = []
        for i in procs_done:
            if i in ['px', 'dn']: suffix_.append(self.fmri_proc_suffix_dict[i])
        suffix = ''.join(suffix_)
        fmri_revPE = self.fmri_revPE_name.replace('.nii', "{}.nii".format(suffix))
        if os.path.isfile(fmri_revPE):
            self.fmri_proc_revPE_name = fmri_revPE
        else:
            assert os.path.isfile(fmri_revPE), ("---'{0}' is nonexistent "
                                                "in {1} !!".format(fmri_revPE, self.fmri_dir) )

        os.chdir(dir_backup)       


            




