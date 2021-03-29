%build functional connectome
function [rho_fc] = piracy_make_connectome(datapath, fmri_name, label_file, output_folder, ifGSR)
    %ifGSR==1: partial correlation with global mean signal as 3rd regressor
    %ifGSR==0: plain pearson's correlation

    ROI_labels = {'ACC L','ACC R','RSC L','RSC R','PPC L','PPC R','MTL L','MTL R','Hip L','Hip R',...
        'Sub L','Sub R','Au L','Au R','V L','V R','S1 L','S1 R','S2 L','S2 R','M L','M R','Str L','Str R',...
        'Tha L','Tha R','HTh L','HTh R'};
    
    nROIs = length(ROI_labels);

    cd(datapath);
    disp(['fMRI file: ', fmri_name]);
    disp(['ROI label file: ', label_file]);
    assert(boolean(exist(label_file, 'file')), label_file);
    if ~exist(output_folder, 'dir'), mkdir(output_folder);end
    %
    scan_ = split(fmri_name, '.nii');
    scan = scan_{1};
    %%
    fnii = load_untouch_nii(fmri_name);
    fmri_data = double(fnii.img);
    nSlices = size(fmri_data, 3);
    brain_mask = (nanmean(fmri_data, 4)>1);
    %
    %% load atlas label2epi nii       
    lnii = load_untouch_nii(label_file);
    label = double(lnii.img);
    ns_label = size(label, 3);  
    assert(nSlices==ns_label, 'num of slices not equal');
    %%
    veclabel = round(label(:));
    label_merge = zeros(size(veclabel));
    %% make merged ROI labels      
    % ACC+CC
    label_merge(ismember(veclabel,[164,165])) = 2; %  R
    label_merge(ismember(veclabel,[212,213])) = 1; %  L
    % RSC
    label_merge(ismember(veclabel,[183,184,185])) = 4; %  R
    label_merge(ismember(veclabel,[231,232,233])) = 3; %  L
    % PPC
    label_merge(ismember(veclabel,[174,178,180,181,182])) = 6; %  R
    label_merge(ismember(veclabel,[222,226,228,229,230])) = 5; %  L
    % MTL
    label_merge(ismember(veclabel,[167,168,170,177,179,197,204])) = 8; %   R
    label_merge(ismember(veclabel,[215,216,218,225,227,245,252])) = 7;  %  L
    % Hip
    label_merge(ismember(veclabel,[73,74,75,76]) ) = 10; %  R
    label_merge(ismember(veclabel,[151,152,153,154])) = 9;  %  L
    % Sub
    label_merge(ismember(veclabel, 70)) = 12;   %  R
    label_merge(ismember(veclabel, 148)) = 11;  %  L
    % AU
    label_merge(ismember(veclabel,[161,162,163])) = 14; %  R
    label_merge(ismember(veclabel,[209,210,211])) = 13; %  L
    % V
    label_merge(ismember(veclabel,[198:203])) = 16; %   R
    label_merge(ismember(veclabel,[246:251])) = 15; %   L
    % S1, S2
    label_merge(ismember(veclabel,[186:195])) = 18; %   R
    label_merge(ismember(veclabel,[234:243])) = 17; %   L
    label_merge(ismember(veclabel, 196)) = 20; %   R
    label_merge(ismember(veclabel, 244)) = 19; %   L
    % M
    label_merge(ismember(veclabel,[175,176])) = 22; % R
    label_merge(ismember(veclabel,[223,224])) = 21; % L

    % CPu striatum
    label_merge(ismember(veclabel,8)) = 24;    %  R 
    label_merge(ismember(veclabel,86)) = 23;   % L 
                
    % Tha       thalamus
    label_merge(ismember(veclabel,17)) = 26;   %  R
    label_merge(ismember(veclabel, 95)) = 25;   %  L
    % HTh       
    label_merge(ismember(veclabel, 26)) = 28;   %  R
    label_merge(ismember(veclabel, 104)) = 27;  %  L  
    % save merged label
    label_merge = reshape(label_merge,size(label));
    nii_label_merge = lnii;
    nii_label_merge.img = int32(label_merge);
    save_untouch_nii(nii_label_merge, replace(label_file, '.nii', '_merge.nii'));    
    %% process fc matrix
    nVols=size(fmri_data, 4);
    %tc_matrix = nan(nVols,nROIs);
    for label_num = 1:nROIs
        label_kb = (label_merge == label_num);
        fmri_roi_k = fmri_data.*label_kb;
        fmri_roi_k(fmri_roi_k == 0) = NaN;
        fmri_roi_m = reshape(fmri_roi_k, [], nVols);
        tc_matrix(label_num,:) = nanmean(fmri_roi_m, 1);
    end
    %
    fmri_data_m = fmri_data.*brain_mask;
    fmri_data_m(fmri_data_m==0)=NaN;
    fmri_data_v = reshape(fmri_data_m, [], nVols);
    m_tc = nanmean(fmri_data_v,1);
    
    %% show time courses of ROI's
    clear mean_tc tc_centered tc_norm 
    mean_tc = nanmean(tc_matrix,2);
    tc_centered = tc_matrix - mean_tc;
    tc_norm = 100 * tc_centered./mean_tc;

    fig1=figure(1);
    set(fig1, 'position', [100, 100, 800, 500]);%set figure size
    subplot(311); imagesc(tc_norm, [-2,2]); colormap(jet)
    axis(gca,'ij');
    set(gca,'YTick',1:nROIs,'YTickLabel',ROI_labels);
    xlabel('normalized ROI timecourses');
    figure(1);
    subplot(312), plot(m_tc);
    xlabel('global average time course');    
    %% functional connectome
    if ifGSR
        %%1. global signal regression
        [rho_fc,pval_fc]=partialcorr(tc_matrix', m_tc', 'Rows', 'complete');%
    else
        [rho_fc,pval_fc]=corrcoef(tc_matrix', 'rows', 'complete');
    end    
   
    z_fc = atanh(rho_fc);
    %%
    figure(1);
    subplot(313); imagesc(rho_fc,[-1 1]); axis square; colormap(jet)
    axis(gca,'ij');
    % Set the remaining axes properties
    set(gca,'Layer','top','XAxisLocation','top','XTick',1:nROIs,'XTickLabel',ROI_labels,'XTickLabelRotation',90);
    scan_name = strrep(scan, '_', '\_');
    suptitle(scan_name)
    set(gca,'YTick',1:nROIs,'YTickLabel',ROI_labels);
    text(15,31,'Functional connectivity','rotation',0,'fontsize',10, ...
        'horizontalalignment','center','verticalalignment','top');
    savefig(fullfile(output_folder, scan));
    saveas(gcf,fullfile(output_folder, [scan, '.jpg']));  
    %%
    figure(2); imagesc(rho_fc,[-1 1]); axis square; colormap(jet)
    axis(gca,'ij');
    % Set the remaining axes properties
    set(gca,'Layer','top','XAxisLocation','top','XTick',1:nROIs,'XTickLabel',ROI_labels,'XTickLabelRotation',90);
    title(scan_name)
    set(gca,'YTick',1:nROIs,'YTickLabel',ROI_labels);

    saveas(gcf,fullfile(output_folder, ['fc_',scan,'.jpg']));
    savefig(fullfile(output_folder, ['fc_',scan]));
    connectome = fullfile(output_folder, ['z_fc_', scan, '.mat']);
    save(connectome,'rho_fc', 'z_fc','pval_fc','tc_matrix', 'ifGSR', 'nSlices', 'nVols', 'fmri_name','label_file', 'ROI_labels'); 

end