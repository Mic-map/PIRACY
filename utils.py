from os import system as scmd
import os
import re
import glob 
import shutil
import nibabel as nib
from tqdm import tqdm
from functools import reduce
import numpy as np
import cv2

HOME_DIR = os.path.dirname(os.path.realpath(__file__))
LIB_DIR = "{}/libs".format(HOME_DIR)

def execute_cmd(cmd, pars_line):
    cmd_line = cmd + " " + pars_line
    scmd(cmd_line)

# create brain mask on fmri image
def fmri_masking(nii_file, suffix='_m', maskSuffix='_mask', createMask=True, applyMask=True, threshold=0.5):
    assert os.path.isfile(nii_file), "File not found!!"
    
    nii_basename = nii_file.split('.nii')[0]
    nii_mean = '{0}_mean'.format(nii_basename)
    mask = '{0}{1}.nii.gz'.format(nii_mean, maskSuffix)   
    nii_tmasked = '' 
    #
    if createMask:        
        execute_cmd('fslmaths', '{0} -Tmean {1}'.format(nii_file, nii_mean))
        # Bias field correction
        execute_cmd('fast', '-t 2 -n 3 -H 0.1 -I 4 -l 20.0 --nopve -B -o {0}'.format(nii_mean))
        execute_cmd('rm', '{0}_seg.nii.gz'.format(nii_mean))

        nii_mean_ = nii_mean
        nii_mean_restore = '{0}_restore.nii.gz'.format(nii_mean)
        if os.path.exists(nii_mean_restore):
            nii_mean_ = nii_mean_restore          
        
        print("create mask based on {0}".format(nii_mean_))        
        execute_cmd('bet', '{0} {1} -f {2} -g 0 -n -m'.format(nii_mean_, nii_mean, threshold))
        # rename
        execute_cmd('mv', '{0}_mask.nii.gz {1}'.format(nii_mean, mask))

    if applyMask:
        suffix_ = "{0}.nii".format(suffix)
        nii_tmasked = nii_file.replace('.nii', suffix_)
        execute_cmd('fslmaths', '{0} -mas {1} {2}'.format(nii_file, mask, nii_tmasked))   
    
    return mask, nii_tmasked


def adjust_mask(mask_name, kernel=(3,3), iterations=1, suffix=''):
    #iterations>0: dilating; <0: eroding; ==0: do nothing
    assert os.path.exists(mask_name), "File not found!!"
    kernel = np.ones(kernel, np.uint8)
    mask_nii = nib.load(mask_name)
    hdr = mask_nii.header
    mask_img = mask_nii.get_fdata()
    nslices = mask_img.shape[2]
    for ns in range(nslices):
        img_ = mask_img[:,:, ns]
        if iterations>0:
            print('dilating mask: ', mask_name)
            mask_img[:,:, ns] = cv2.dilate(img_, kernel, iterations)
        elif iterations<0:
            print('eroding mask: ', mask_name)
            mask_img[:,:, ns] = cv2.erode(img_, kernel, int(iterations))
        else:
            return

    mask_data = nib.Nifti1Image(mask_img, None, header=hdr)
    if suffix: mask_name = mask_name.replace('.nii', '{0}.nii'.format(suffix))
    nib.save(mask_data, mask_name) 

def fmri_topup(datapath, fmri_name, fmri_revPE, suffix='_uw', tpFolder='Topup', acqp='acq_pars.txt', mycnf='mycnf_x10.cnf'):
    print('Start topup correction for {0}...'.format(fmri_name))
    orgDir = os.getcwd()
    os.chdir(datapath)
    if not os.path.exists(tpFolder): os.mkdir(tpFolder)
    assert os.path.isfile(fmri_name), "File not found!!"
    assert os.path.isfile(fmri_revPE), "File not found"

    fmri_nii = nib.load(fmri_name)
    fmri_img = fmri_nii.get_fdata()
    linPE_hdr = fmri_nii.header
    dims = linPE_hdr['dim']
    slices = dims[3]#number of slices
    print('Number of slices: ', slices)

    fniiName = fmri_name.split('.nii')[0] #base name without extension
    fniiName_uw = "{0}{1}".format(fniiName, suffix)
    fmri_uw = fniiName_uw + '.nii.gz'

    #enter Topup folder {
    os.chdir(tpFolder)
    #create a temporary short linPE fmri run
    linPE_img = fmri_img[:,:,:,:10]
    linPE_hdr['dim'][4]=10
    linPE_nii = nib.nifti1.Nifti1Image(linPE_img, None, header=linPE_hdr)
    fmri_linPE = "tmp_fmri_short_linPE.nii.gz"
    nib.save(linPE_nii, fmri_linPE)
    
    execute_cmd("fslroi","{0} epi_lin_run 8 2".format(fmri_linPE))

    execute_cmd("fslroi","../{0} epi_rev_run 8 2".format(fmri_revPE))

    execute_cmd("fslmerge","-t epi_tu_run epi_lin_run.nii.gz epi_rev_run.nii.gz")

    execute_cmd("fslroi","epi_tu_run.nii.gz epi_tu_run_first 0 -1 0 -1 0 1 0 -1")

    execute_cmd("fslroi","epi_tu_run.nii.gz epi_tu_run_last 0 -1 0 -1 {0} 1 0 -1".format(slices-1))

    execute_cmd("fslmerge", "-z epi_tu_run_padded epi_tu_run_first.nii.gz" \
        " epi_tu_run.nii.gz epi_tu_run_last.nii.gz")
    
    print("Parameter file: {}".format(mycnf))

    execute_cmd("topup","--imain=epi_tu_run_padded.nii.gz" \
    " --datain={0}/{2}" \
    " --out=topup_results_run" \
    " --iout=epi_tu_run_padded_uw" \
    " --config={0}/{1} -v".format(LIB_DIR, mycnf, acqp)) 

    execute_cmd("fslroi", "../{0} {1}_first 0 -1 0 -1 0 1 0 -1".format(fmri_name, fniiName))

    execute_cmd("fslroi", "../{0} {1}_last 0 -1 0 -1 {2} 1 0 -1".format(fmri_name, fniiName, slices-1))

    execute_cmd("fslmerge", "-z {0}_padded {0}_first.nii.gz" \
        " ../{1} {0}_last.nii.gz".format(fniiName, fmri_name))

    execute_cmd("applytopup", "-i {0}_padded.nii.gz --inindex={1}" \
            " --datain={2}/{3}" \
            " --topup=topup_results_run --out={0}_padded_uw" \
            " --method=jac".format(fniiName, 1, LIB_DIR, acqp))

    execute_cmd("fslroi", "{0}_padded_uw.nii.gz" \
        " ../{1} 0 -1 0 -1 1 {2} 0 -1".format(fniiName, fniiName_uw, slices))

    assert os.path.exists('../{0}'.format(fmri_uw)), "Error----Topup failed!!"

    execute_cmd("rm", "epi_lin_run.nii.gz epi_rev_run.nii.gz" \
        " epi_tu_run_first.nii.gz epi_tu_run_last.nii.gz" \
        " epi_tu_run_padded.nii.gz" \
        " {0}_first.nii.gz {0}_last.nii.gz" \
        " {0}_padded_uw.nii.gz".format(fniiName))

    os.chdir(orgDir)

    return fmri_uw

def ants_reg_anat_to_atlas(anat_dir, anat_name, atlasLabel='nmd2WHS_label_v8_ds_ForfMRI_x10.nii', atlasTemp='nmd2WHS_brain_v2_ds_x10.nii.gz',
                                regFolder='reg', atlasDir='', output_suffix='_anat2tempA_', label_suffix='_label'):
    print('-----Registering {0} to {1}'.format(anat_name, atlasTemp))
    orgDir = os.getcwd()
    os.chdir(anat_dir)

    if not os.path.isdir(regFolder): os.mkdir(regFolder)

    if atlasDir:
        atlas_dir = atlasDir
    else:
        atlas_dir= LIB_DIR
    
    os.chdir(regFolder) #enter reg folder
    moving_nii = "../{}".format(anat_name)
    fix_nii='{0}/{1}'.format(atlas_dir,atlasTemp)

    if atlasLabel:
        label='{0}/{1}'.format(atlas_dir, atlasLabel)
        assert os.path.exists(label), "---Atlas label not found!!"

    anat_basename = anat_name.split('.nii')[0]
    OUTPUT_PREFIX='{0}{1}'.format(anat_basename, output_suffix)

    execute_cmd("antsRegistration",  "--dimensionality 3 --float 0 --interpolation BSpline" \
		         " --use-histogram-matching 1" \
		         " --winsorize-image-intensities [0.005,0.995]" \
		         " --output [{0},{0}Warped.nii.gz,{0}InverseWarped.nii.gz]" \
		         " --initial-moving-transform [{1},{2},1]" \
                 " --transform translation[0.1]" \
		         "   --metric MI[{1},{2},1,32,Random,0.25]" \
		         "   --convergence [50,1e-6,10]" \
		         "   --shrink-factors 1" \
		         "   --smoothing-sigmas 0vox" \
		         " --transform Rigid[0.05]" \
		         "   --metric MI[{1},{2},1,64,Random,0.25]" \
		         "   --convergence [500x250x50,1e-6,10]" \
		         "   --shrink-factors 2x2x1" \
		         "   --smoothing-sigmas 2x1x0vox" \
		         " --transform SyN[0.1,3,0]" \
		         "   --metric CC[{1},{2},1,4]" \
		         "   --convergence [50x10,1e-6,10]" \
		         "   --shrink-factors 2x1" \
		         "   --smoothing-sigmas 1x0vox".format(OUTPUT_PREFIX, fix_nii, moving_nii))    

    #" --input-image-type 2" \
    execute_cmd("antsApplyTransforms", \
            "--dimensionality 3" \
            " --input {0}"\
            " --reference-image {1}" \
            " --transform [{2}0GenericAffine.mat, 1]" \
            " --transform {2}1InverseWarp.nii.gz" \
            " --output ../{3}{4}.nii.gz" \
            " --interpolation NearestNeighbor".format(label, moving_nii, OUTPUT_PREFIX, anat_basename, label_suffix)) 

    os.chdir(orgDir)
    return "{0}{1}.nii.gz".format(anat_basename, label_suffix)

def ants_reg_fmri_to_anat(fmri_dir, anat_brain_file, fmri_name, anatLabel, label_suffix='_label', regFolder='reg', applyOnly=False):	
	orgDir = os.getcwd()
	os.chdir(fmri_dir)
	if not os.path.exists(regFolder): os.mkdir(regFolder)
	assert os.path.isfile(anatLabel), "----Anatomical label file not found!!"
	
	os.chdir(regFolder)
	
	fmriMean = fmri_name.replace('.nii', '_mean.nii')
	print('Creating {0} ...'.format(fmriMean))
	execute_cmd('fslmaths', '../{0} -Tmean {1}'.format(fmri_name, fmriMean))
	
	moving_nii = fmriMean
	fix_nii = anat_brain_file
	fix_label = anatLabel
	
	fmri_bsname = fmri_name.split('.nii')[0]
	output_prefix='{0}_2anat_'.format(fmri_bsname)
			
	if not applyOnly:
		###
		print("reg {0} to {1}".format(moving_nii, fix_nii))
		execute_cmd("antsRegistration", " --dimensionality 3 --float 0 --interpolation BSpline" \
				" --use-histogram-matching 1" \
				" --winsorize-image-intensities [0.005,0.995]" \
				" --output [{0},{0}Warped.nii.gz,{0}InverseWarped.nii.gz]" \
				" --initial-moving-transform [{1},{2},1]" \
				" --transform translation[0.1]" \
				"  --metric MI[{1},{2},1,32,Random,0.25]" \
				"  --convergence [50,1e-6,10]" \
				"  --shrink-factors 1" \
				"  --smoothing-sigmas 0vox" \
				" --transform Rigid[0.05]" \
				"  --metric MI[{1},{2},1,64,Random,0.25]" \
				"  --convergence [500x250x50,1e-6,10]" \
				"  --shrink-factors 2x2x1" \
				"  --smoothing-sigmas 2x1x0vox" \
				" --transform SyN[0.1,3,0]" \
				"  --metric CC[{1},{2},1,6]" \
				"  --convergence [50x10,1e-6,10]" \
				"  --shrink-factors 2x1" \
				"  --smoothing-sigmas 1x0vox".format(output_prefix, fix_nii, moving_nii))
	
	#map anat label to the fmri space
	fmri_label = fmri_name.replace('.nii', '{}.nii'.format(label_suffix))
	print("applying warp")#
	execute_cmd("antsApplyTransforms", " --dimensionality 3" \
			" --input {0}" \
			" --reference-image {1}" \
			" --transform [{2}0GenericAffine.mat, 1]" \
			" --transform {2}1InverseWarp.nii.gz" \
			" --output ../{3}" \
			" --interpolation NearestNeighbor".format(fix_label, moving_nii, output_prefix, fmri_label)) 
	os.chdir(orgDir)  
	
	return fmri_label

def fmri_melodic_ica(datapath, fmri_name, outComponents=0, useFix=True, cutoff=100, ndelete=10, mcflirt=0, template='template_fsems_anat_restore_brain.nii'):
    orgDir = os.getcwd()
    #outComponent <= 0: enable a utomatic dimensionality estimation
    auto_dim = 1 if outComponents <= 1 else 0
    if useFix:
        tptPath = os.path.join(LIB_DIR, template)
        assert os.path.exists(tptPath)
        tptPathStrs = '\/'.join(LIB_DIR.split('/') )

    os.chdir(datapath)
    PathNii = os.path.abspath(fmri_name)
    niiStr = '\/'.join(PathNii.split('/'))

    info = nib.load(fmri_name).header
    tr = float(info['pixdim'][4])
    totalVoxels = 1
    dims=info['dim']
    rp = int(dims[4])
    for d in dims[1:5]:
        totalVoxels = totalVoxels*d

    if useFix:
        execute_cmd("sed", " \"s/\\/data\\/saepi_fmri_05_dn_uw_m/" \
                "{0}/g\"" \
                " {1}/feat_design_fix.fsf" \
                " > feat_design.fsf".format(niiStr,LIB_DIR))
        execute_cmd("sed", " -i \"s/\\/template\\/templateBrain/{0}\\/{1}/g\" feat_design.fsf".format(tptPathStrs,template))
    else:
        execute_cmd("sed", " \"s/\\/data\\/saepi_fmri_05_dn_uw_m/" \
                "{0}/g\"" \
                " {1}/feat_design_no_fix.fsf" \
                " > feat_design.fsf".format(niiStr,LIB_DIR))

    execute_cmd("sed", " -i \"s/fmri(tr) 1.500000/fmri(tr) {0}/g\" feat_design.fsf".format(tr))
    execute_cmd("sed", " -i \"s/fmri(npts) 400/fmri(npts) {0}/g\" feat_design.fsf".format(rp))
    execute_cmd("sed", " -i \"s/fmri(totalVoxels) 14745600/fmri(totalVoxels) {0}/g\" feat_design.fsf".format(totalVoxels))
    # High pass filter cutoff
    execute_cmd("sed", " -i \"s/fmri(paradigm_hp) 100/fmri(paradigm_hp) {0}/g\" feat_design.fsf".format(cutoff))
    #number of volumes to delete
    execute_cmd("sed", " -i \"s/fmri(ndelete) 10/fmri(ndelete) {0}/g\" feat_design.fsf".format(ndelete))

    #set output components
    if auto_dim==0:
        execute_cmd("sed", " -i \"s/fmri(dim) 40/fmri(dim) {0}/g\" feat_design.fsf".format(outComponents))
        print("number of components: ", outComponents)

    #Automatic dimensionality estimation
    execute_cmd("sed", " -i \"s/fmri(dim_yn) 0/fmri(dim_yn) {0}/g\" feat_design.fsf".format(auto_dim))
    #motion correction or not
    execute_cmd("sed", " -i \"s/fmri(mc) 0/fmri(mc) {0}/g\" feat_design.fsf".format(mcflirt))
    execute_cmd("feat","feat_design.fsf")

    os.chdir(orgDir)

def create_dummy_mcf_par(fmri_ica):
    #generate dummy motion parameter file for FIX
    assert os.path.isdir(fmri_ica), "---ica folder not found!!"
    orgDir = os.getcwd()
    os.chdir(fmri_ica)

    assert os.path.exists("filtered_func_data.nii.gz"), "----'filtered_func_data.nii.gz' not found!!"
    info = nib.load("filtered_func_data.nii.gz").header
    nrep=info['dim'][4]
    if not os.path.exists('mc'): os.mkdir('mc')
    with open("mc/prefiltered_func_data_mcf.par", 'w') as f:
        f.write('-1.81077e-08  -0  1.18916e-08  0  0  1.3681e-08 \n')
        for n in range(nrep-2):
            f.write('0  -0  0  0  0  0   \n')
        f.write('-1.81077e-08  -0  1.18916e-08  0  0  1.3681e-08 ')

    os.chdir(orgDir)
    return True

def fix_ica_classify(datapath, ica, threshold, model=''):
    #classify noise ICA components
    if not model:
        model = '{0}/Fix_base_ic40_v2.RData'.format(LIB_DIR)
    assert os.path.exists(model), "---FIX model file not found!!"
    modelName = os.path.basename(model)
    modelName = os.path.splitext(modelName)[0]	        
    
    orgDir = os.getcwd()
    os.chdir(datapath)   
    assert os.path.exists(ica), "---ICA folder not found!!"  
                                                                         
    resultFile = "fix4melview_{0}_thr{1}.txt".format(modelName, threshold)
    resultFile = os.path.join(ica, resultFile)

    if not os.path.exists(resultFile):
        print('fix -c {0} {1} {2}'.format(ica, model, threshold))
        execute_cmd('fix', '-c {0} {1} {2}'.format(ica, model, threshold))
    
    os.chdir(orgDir)
    resultFile = os.path.join(datapath, resultFile)

    return resultFile

def fmri_ica_clean_noise(datapath, ica, noise_file='hand_label_noise.txt', suffix='_clean'): 
    orgDir = os.getcwd()
    os.chdir(datapath)                                 

    assert os.path.exists(ica), ica
    noise_file = os.path.join(ica, noise_file)
    assert os.path.exists(noise_file), "---Noise ICA file not found:{}!!".format(noise_file)
    noise_list = read_noiseList_file(noise_file, splitter=' ')
    ica_bsname = os.path.splitext(ica)[0]
    cleaned_bsname = ''.join([ica_bsname, suffix])
    cleaned = "{0}.nii.gz".format(cleaned_bsname)

    if noise_list:
       os.chdir(ica)
       noise_str = ','.join(map(lambda i:str(i), noise_list)) #convert to string: "1, 2, 3 ..."
       print('noise list: ')
       print(noise_str)
       execute_cmd("fsl_regfilt", " -i filtered_func_data.nii.gz -o ../{0}" \
                                  " -d filtered_func_data.ica/melodic_mix -f \"{1}\"".format(cleaned_bsname, noise_str))
       print('Cleaned data: {0}'.format(cleaned))
    else:
       print("No noise list!!!")

    os.chdir(orgDir)

    return cleaned

def read_noiseList_file(predFile, splitter=','):
    rexp = re.compile('\[([\S\s]*)\]')
    with open(predFile, 'r') as f:
        content = f.read()
    pred = re.findall(rexp, content)[0]
    if not pred:
        print("no noise component found!")
        return []

    pred = pred.split(splitter)
    if splitter == ',':
        predList = [int(i) for i in pred if i.strip()]
    else:
        predList=[]
        for i in pred:
            if i.strip():
                predList.append(int(i.strip().split(',')[0]))
    #print("Noise list: ", predList)
    return predList