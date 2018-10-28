% do QSM (RESHARP + inversion) on individual echoes

% (1) test bestpath unwrapping
mkdir('IND/bestpath')
cd IND/bestpath
R = 1;
for i = 1:imsize(4)
	tfs(:,:,:,i) = unph(:,:,:,i)/TE(i);
end

if (strcmp(readout,'unipolar'))
    tfs = tfs/(2.675e8*dicom_info.MagneticFieldStrength)*1e6; % unit ppm
else 
    tfs = -tfs/(2.675e8*dicom_info.MagneticFieldStrength)*1e6; % unit ppm
end

nii = make_nii(tfs,vox);
save_nii(nii,'tfs.nii');



for i = 1:imsize(4)
	[lfs_resharp(:,:,:,i), mask_resharp] = resharp(tfs(:,:,:,i),mask.*R,vox,smv_rad,tik_reg,cgs_num);
    % % 3D 2nd order polyfit to remove any residual background
    % lfs_resharp= lfs_resharp - poly3d(lfs_resharp,mask_resharp);

    % save nifti
    [~,~,~] = mkdir('RESHARP');
    nii = make_nii(lfs_resharp(:,:,:,i),vox);
    save_nii(nii,['RESHARP/e' num2str(i) '_lfs_resharp_tik_', num2str(tik_reg), '_num_', num2str(cgs_num), '.nii']);

    % % inversion of susceptibility 
    % disp('--> TV susceptibility inversion on RESHARP...');
    % sus_resharp = tvdi(lfs_resharp,mask_resharp,vox,tv_reg,mag(:,:,:,end),z_prjs,inv_num); 
   
    % % save nifti
    % nii = make_nii(sus_resharp.*mask_resharp,vox);
    % save_nii(nii,['RESHARP/sus_resharp_tik_', num2str(tik_reg), '_tv_', num2str(tv_reg), '_num_', num2str(inv_num), '.nii']);
    
    
    % iLSQR
    chi_iLSQR(:,:,:,i) = QSM_iLSQR(lfs_resharp(:,:,:,i)*(2.675e8*dicom_info.MagneticFieldStrength)/1e6,mask_resharp,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
    nii = make_nii(chi_iLSQR(:,:,:,i),vox);
    save_nii(nii,['RESHARP/e' num2str(i) '_chi_iLSQR_smvrad' num2str(smv_rad) '.nii']);
    
 %    % MEDI
 %    %%%%% normalize signal intensity by noise to get SNR %%%
 %    %%%% Generate the Magnitude image %%%%
 %    if imsize(4) > 1
 %        iMag = sqrt(sum(mag.^2,4));
	% else
	% 	iMag = mag;
 %    end
    
 %    % [iFreq_raw N_std] = Fit_ppm_complex(ph_corr);
 %    matrix_size = single(imsize(1:3));
 %    voxel_size = vox;
 %    if imsize(4) > 1
 %        delta_TE = TE(2) - TE(1);
 %    else
 %        delta_TE = 0.001;
 %    end
 %    B0_dir = z_prjs';
 %    CF = dicom_info.ImagingFrequency *1e6;
 %    iFreq = [];
 %    N_std = 1;
 %    RDF = lfs_resharp(:,:,:,i)*2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6;
 %    Mask = mask_resharp;
 %    save RDF.mat RDF iFreq iMag N_std Mask matrix_size...
 %         voxel_size delta_TE CF B0_dir;
 %    QSM = MEDI_L1('lambda',1000);
 %    nii = make_nii(QSM.*Mask,vox);
 %    save_nii(nii,['RESHARP/e' num2str(i) '_MEDI1000_RESHARP_smvrad' num2str(smv_rad) '.nii']);
 %    QSM = MEDI_L1('lambda',2000);
 %    nii = make_nii(QSM.*Mask,vox);
 %    save_nii(nii,['RESHARP/e' num2str(i) '_MEDI1000_RESHARP_smvrad' num2str(smv_rad) '.nii']);
end





% (2) test laplacian unwrapping and QSM of each echo
mkdir('IND/laplacian')
cd IND/laplacian

    disp('--> unwrap aliasing phase using laplacian...');
    Options.voxelSize = vox;
    for i = 1:imsize(4)
        unph(:,:,:,i) = lapunwrap(ph_corr(:,:,:,i), Options).*mask;
    end

    nii = make_nii(unph, vox);
    save_nii(nii,'unph_lap.nii');

    unph_lap = unph;
    % save('raw.mat','unph_lap','-append');




for i = 1:imsize(4)
	tfs(:,:,:,i) = unph(:,:,:,i)/TE(i);
end

if (strcmp(readout,'unipolar'))
    tfs = tfs/(2.675e8*dicom_info.MagneticFieldStrength)*1e6; % unit ppm
else 
    tfs = -tfs/(2.675e8*dicom_info.MagneticFieldStrength)*1e6; % unit ppm
end

nii = make_nii(tfs,vox);
save_nii(nii,'tfs.nii');



R = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



nii = load_nii('e1_lfs_resharp_tik_0.0001_num_200.nii');
e1_lfs = single(nii.img);
nii = load_nii('e2_lfs_resharp_tik_0.0001_num_200.nii');
e2_lfs = single(nii.img);
nii = load_nii('e3_lfs_resharp_tik_0.0001_num_200.nii');
e3_lfs = single(nii.img);
nii = load_nii('e4_lfs_resharp_tik_0.0001_num_200.nii');
e4_lfs = single(nii.img);

lfs_ave3 = (e2_lfs + e3_lfs + e4_lfs)/3;
nii = make_nii(lfs_ave3);
save_nii(nii,'lfs_ave3.nii');

lfs_ave4 = (e1_lfs + e2_lfs + e3_lfs + e4_lfs)/4;
nii = make_nii(lfs_ave4);
save_nii(nii,'lfs_ave4.nii');



nii = load_nii('e1_chi_iLSQR_smvrad3.nii');
e1_chi = single(nii.img);
nii = load_nii('e2_chi_iLSQR_smvrad3.nii');
e2_chi = single(nii.img);
nii = load_nii('e3_chi_iLSQR_smvrad3.nii');
e3_chi = single(nii.img);
nii = load_nii('e4_chi_iLSQR_smvrad3.nii');
e4_chi = single(nii.img);


chi_ave3 = (e2_chi + e3_chi + e4_chi)/3;
nii = make_nii(chi_ave3);
save_nii(nii,'chi_ave3.nii');

chi_ave4 = (e1_chi + e2_chi + e3_chi + e4_chi)/4;
nii = make_nii(chi_ave4);
save_nii(nii,'chi_ave4.nii');













