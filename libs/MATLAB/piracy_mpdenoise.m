%% by Yujian Diao, Denoise for fMRI data
% run on linux server, called by python
% Modified script from Ileana Jelescu 'rs_fmri_MPdenoising.m'
% 'MPdenoising.m by Jelle Veraart
% https://github.com/NYU-DiffusionMRI/mppca_denoise/blob/master/MPdenoising.m

%%
function ret_niis = piracy_mpdenoise(dataDir, fmri_nii, fmri_revPE_nii, fmri_mask, suffix, saveAll)
    %
    ret_niis = {};
    cd(dataDir);

    %%
    kernel = [5 5 5];
    nrep_list=[]; 
    epi_img_list = {};
    epi_hdr_list = {};

    fMRI_epis = {fmri_nii, fmri_revPE_nii};
    for fn=1: length(fMRI_epis)
        fname=fMRI_epis{fn};
        disp(fname);

        epi_nii = load_untouch_nii(fname); 
        epi_hdr_list{fn} = epi_nii.hdr;
        dims = epi_nii.hdr.dime;
        epi_img_list{fn} = epi_nii.img; 
        Nimg = size(epi_nii.img, 4); clear epi_nii;  
        nrep_list(fn)=Nimg;

        if fn==1           
            ro_pixdim = dims.pixdim(2);
            pe_pixdim = dims.pixdim(3); 
            thk_epi = dims.pixdim(4);
        end

    end
    disp(nrep_list);

    %%
    vol_2 = [];
    for ne=1:length(fMRI_epis)
        vol_2 = cat(4,vol_2,epi_img_list{ne});    
    end

    disp(['Using mask: ', fmri_mask]);
    mask1 = load_untouch_nii(fmri_mask); 
    hdr_mask = mask1.hdr;
    mask = double(mask1.img);    
    %slice padding
    buffer_m = cat(3,zeros(size(mask,1),size(mask,2),2),mask,zeros(size(mask,1),size(mask,2),3));
    buffer_m = buffer_m>0;
    buffer = cat(3,vol_2(:,:,end-1:end,:),vol_2,vol_2(:,:,1:3,:));

    disp('Denoising...');
    tic;
    [signal,sigma, P]=MPdenoising(buffer,buffer_m,kernel,'full');
    toc 
    signal_unpadded = signal(:,:,3:end-3,:);
    sigma_unpadded = sigma(:,:,3:end-3,:);
    P_unpadded = P(:,:,3:end-3,:);
 
    for ifn=1:length(fMRI_epis)
        signalx = signal_unpadded(:,:,:,1:nrep_list(ifn));
        nii = make_nii(signalx,[ro_pixdim pe_pixdim thk_epi],[0 0 0],16);
        nii.hdr = epi_hdr_list{ifn};
        %nii.hdr.dime.pixdim(5) = tr;
        newName = replace(fMRI_epis{ifn}, '.nii', [suffix, '.nii']);
        save_nii(nii,newName)
        ret_niis{ifn} = newName;
        if saveAll
            nii = make_nii(signalx./sigma_unpadded,[ro_pixdim pe_pixdim thk_epi],[0 0 0],16);
            nii.hdr = epi_hdr_list{ifn};
            save_nii(nii,replace(fMRI_epis{ifn}, '.nii', '_snr.nii'))
        end
        if ~isempty(signal_unpadded)
            signal_unpadded(:,:,:,1:nrep_list(ifn))=[];
        end
    end 

    if saveAll
        nii = make_nii(sigma_unpadded,[ro_pixdim pe_pixdim thk_epi],[0 0 0],16);
        if ~isempty(fmri_mask), nii.hdr = hdr_mask; end 
        nii.hdr.dime.datatype=16;
        save_nii(nii,replace(fMRI_epis{1}, '.nii', '_sigma.nii'));

        nii = make_nii(P_unpadded,[ro_pixdim pe_pixdim thk_epi],[0 0 0],16);
        if ~isempty(fmri_mask), nii.hdr = hdr_mask; end
        nii.hdr.dime.datatype=16;
        save_nii(nii,replace(fMRI_epis{1}, '.nii', '_P.nii'));
    end

end
    

