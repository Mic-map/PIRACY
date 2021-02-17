%scaling pixdim x 10 in nii_header 
function scale_file = piracy_nii_scale_pixdim_x10(nii_directory, nii_name, suffix)
    %suffix= '_x10.nii'
    disp(['Scale pixdim by 10 in the header of ', nii_name, '']);
    cd(nii_directory);
    %% 1) load data
    nii = load_untouch_nii(nii_name);
    nii_img = double(nii.img);
    tr = nii.hdr.dime.pixdim(5);
      
    scale_file = replace(nii_name,'.nii', suffix);
    disp(['Generate ', scale_file]);
    nii2 = make_nii(nii_img,[nii.hdr.dime.pixdim(2:4)*10],[0 0 0],16); 
    nii2.hdr.dime.pixdim(5) = tr;    % TR
    save_nii(nii2,scale_file);
  
end
