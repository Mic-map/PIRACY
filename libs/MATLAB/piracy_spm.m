%% Slice time correction on topup corrected, brain extracted images
% 1) Load nifti fMRI data (xxx_dn_mc.nii or xxx_MP)
% 2) Slice timing correction: spm_slice_timing.m
%% changes, smoothing kernel = 0.5 in throgh plane direction
% TR is always 1.6seconds -- 0.8s means two shots EPI
% TA -- acquisition time, TA = TR - TR/num_Slice
% sliceOrder -- Interleaved Ascending [ 1 3 5 7 9 2 4 6 8]
% Reference slice -- 1st slice
% timing: [time between slices, time between last slice and next volume]
% smooth.fwhm = [0.36 0.36 1/2]*10 -- smooth kernel, no smoothing through slices
% smooth.im = 0 -- no Implicit masking
% prefix: output file starts with 'a' for STC, 'sa' for smoothing
%% access path with data
function ret_name = piracy_spm(dataDir,  file_name_gz, suffix)%prefix_stc, prefix_sm
    cd(dataDir);

    disp(['Pre-process for (', file_name_gz, ')']);
    if contains(file_name_gz, '.gz')
        file_name = erase(file_name_gz, '.gz');
            disp(['gunzipping file: ', file_name_gz]);
            gunzip(file_name_gz);
    else
        file_name = file_name_gz;
    end    

    %1) load data, create binary mask based on current data
    nii_epi_4D = load_untouch_nii(file_name);
    TR = nii_epi_4D.hdr.dime.pixdim(5);
    sm_kernel = nii_epi_4D.hdr.dime.pixdim(2:4);
    epi_4D = double(nii_epi_4D.img);
    frame = size(epi_4D,4); num_Slice = size(epi_4D,3);
    mask = epi_4D(:,:,:,1)>1;  
    
    %2) Slice timing correction
    file = file_name;
    % Interleaved Ascending
    sliceOrder = [1:2:num_Slice 2:2:num_Slice]; 
    disp("slice order:");
    disp(sliceOrder);
    
    refSlice = 1;
    TA = TR - TR/num_Slice;
    timing = [TA/(num_Slice-1), TR-TA];
    %
    file_STC = ['a',file];
    %file_STC = replace(file, '.nii', '_a.nii');
    stc_gz = replace(file_STC, '.nii', '.nii.gz');
    if exist(stc_gz, 'file')        
        delete(stc_gz);
    end    
    spm_slice_timing(file, sliceOrder, refSlice, timing, 'a');
    
    %3) Spatial smoothing
    disp('Spatial smoothing...');
    %smooth_filename = [prefix_sm, prefix_stc, file_name];
    smooth_filename = replace(file_name, '.nii', [suffix, '.nii']);
    dtype = 0; % same data type
    spm_smooth(file_STC, smooth_filename, sm_kernel, dtype);
    
    %4) maskout boundary pixels
    epi_4D_seg = zeros(size(epi_4D));
    sa_epi = load_untouch_nii(smooth_filename);
    sa_epi_4D = double(sa_epi.img);
    for j = 1 : frame
        for ss = 1: num_Slice
            epi_4D_seg(:,:,ss,j) = sa_epi_4D(:,:,ss,j) .* mask(:,:,ss);
        end
    end

    nii = nii_epi_4D;
    nii.img = epi_4D_seg;
    ret_name = replace(smooth_filename, '.nii', '.nii.gz');
    save_untouch_nii(nii,ret_name);
    
    %5)
    if exist(file_name_gz, 'file')
        delete(erase(file_name_gz, '.gz'));
    end

    %gzip(file_STC);
    delete(file_STC);
end
