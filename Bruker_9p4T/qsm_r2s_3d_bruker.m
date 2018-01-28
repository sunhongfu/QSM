function qsm_r2s_3d_bruker(path_fid,path_out)

[iField voxel_size matrix_size TE delta_TE CF Affine3D B0_dir TR NumEcho] = Read_Bruker_raw_sun(path_fid);

% define directories
path_qsm = [path_out '/QSM_R2s_9p4'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);

% need to be large so that can use N3 correction
% how about change to N4 (ants?)
nii = make_nii(abs(iField)*10000,voxel_size); 
save_nii(nii,'mag_all.nii');
nii = make_nii(angle(iField),voxel_size);
save_nii(nii,'ph_all.nii');

imsize = size(iField);
z_prjs = B0_dir;

mkdir src
for i = 1:imsize(4)
	nii = make_nii(abs(iField(:,:,:,i))*10000,voxel_size);
	save_nii(nii,['src/mag' num2str(i) '.nii']);

	% N3 correction
	setenv('echonum',num2str(i));
	unix('nii2mnc src/mag${echonum}.nii src/mag${echonum}.mnc');
	unix('nu_correct src/mag${echonum}.mnc src/corr_mag${echonum}.mnc -V1.0 -distance 10');
	unix('mnc2nii src/corr_mag${echonum}.mnc src/corr_mag${echonum}.nii');

end
!rm src/*.mnc
!rm src/*.imp

ph = angle(iField);

for i = 1:imsize(4)
	nii = make_nii(ph(:,:,:,i),voxel_size);
	save_nii(nii,['src/ph' num2str(i) '.nii']);
end

mag = zeros(imsize);
for echo = 1:imsize(4)
    nii = load_nii(['src/corr_mag' num2str(echo) '.nii']);
    mag(:,:,:,echo) = double(nii.img);
end

save raw.mat


% % 1 manual masking
% % extract brain using itk-snap
% nii = load_nii('mask.nii');
% mask = double(nii.img);

% 2 thresholding magnitude of TE1
% mask based on mag1
mask = ones(size(mag(:,:,:,1)));
mag1_blur = smooth3(mag(:,:,:,1),'box',round(0.1./voxel_size)*4+1); 
mask(mag1_blur<4000) = 0;
nii = make_nii(mask,voxel_size);
save_nii(nii,'mask_thr.nii');



% phase offset correction
ph_corr = geme_cmb_mouse(mag.*exp(1j*ph),voxel_size,TE,mask);
% save offset corrected phase niftis
for echo = 1:imsize(4)
    nii = make_nii(ph_corr(:,:,:,echo),voxel_size);
    save_nii(nii,['src/corr_ph' num2str(echo) '.nii']);
end


% % unwrap each echo using prelude (too slow)
% disp('--> unwrap aliasing phase for all TEs using prelude...');
% setenv('echo_num',num2str(size(iField,4)));
% bash_command = sprintf(['for ph in src/corr_ph[1-$echo_num].nii\n' ...
% 'do\n' ...
% '   base=`basename $ph`;\n' ...
% '   dir=`dirname $ph`;\n' ...
% '   mag=$dir/"corr_mag"${base:7};\n' ...
% '   unph="unph"${base:7};\n' ...
% '   prelude -a $mag -p $ph -u $unph -m mask_thr.nii -n 12&\n' ...
% 'done\n' ...
% 'wait\n' ...
% 'gunzip -f unph*.gz\n']);
% unix(bash_command);


% unph = zeros(imsize);
% for echo = 1:imsize(4)
%     nii = load_nii(['unph' num2str(echo) '.nii']);
%     unph(:,:,:,echo) = double(nii.img);
% end



%%%%%%% bestpath unwrapping (fast)
% unwrap the phase using best path
disp('--> unwrap aliasing phase using bestpath...');
mask_unwrp = uint8(abs(mask)*255);
fid = fopen('mask_unwrp.dat','w');
fwrite(fid,mask_unwrp,'uchar');
fclose(fid);

[pathstr, ~, ~] = fileparts(which('3DSRNCP.m'));
setenv('pathstr',pathstr);
setenv('nv',num2str(imsize(1)));
setenv('np',num2str(imsize(2)));
setenv('ns',num2str(imsize(3)));

unph = zeros(imsize);

for echo_num = 1:imsize(4)
    setenv('echo_num',num2str(echo_num));
    fid = fopen(['wrapped_phase' num2str(echo_num) '.dat'],'w');
    fwrite(fid,ph_corr(:,:,:,echo_num),'float');
    fclose(fid);
    if isdeployed
        bash_script = ['~/bin/3DSRNCP wrapped_phase${echo_num}.dat mask_unwrp.dat ' ...
        'unwrapped_phase${echo_num}.dat $nv $np $ns reliability${echo_num}.dat'];
    else    
        bash_script = ['${pathstr}/3DSRNCP wrapped_phase${echo_num}.dat mask_unwrp.dat ' ...
        'unwrapped_phase${echo_num}.dat $nv $np $ns reliability${echo_num}.dat'];
    end
    unix(bash_script) ;

    fid = fopen(['unwrapped_phase' num2str(echo_num) '.dat'],'r');
    tmp = fread(fid,'float');
    % tmp = tmp - tmp(1);
    unph(:,:,:,echo_num) = reshape(tmp - round(mean(tmp(mask==1))/(2*pi))*2*pi ,imsize(1:3)).*mask;
    fclose(fid);

    fid = fopen(['reliability' num2str(echo_num) '.dat'],'r');
    reliability_raw = fread(fid,'float');
    reliability_raw = reshape(reliability_raw,imsize(1:3));
    fclose(fid);

    nii = make_nii(reliability_raw.*mask,voxel_size);
    save_nii(nii,['reliability_raw' num2str(echo_num) '.nii']);
end

nii = make_nii(unph,voxel_size);
save_nii(nii,'unph_bestpath.nii');



% check and correct for 2pi jump between echoes
disp('--> correct for potential 2pi jumps between TEs ...')


nii = load_nii('unph_diff.nii');
unph_diff = double(nii.img);

for echo = 2:imsize(4)
    meandiff = unph(:,:,:,echo)-unph(:,:,:,1)-double(echo-1)*unph_diff;
    meandiff = meandiff(mask==1);
    meandiff = mean(meandiff(:));
    njump = round(meandiff/(2*pi));
    disp(['    ' num2str(njump) ' 2pi jumps for TE' num2str(echo)]);
    unph(:,:,:,echo) = unph(:,:,:,echo) - njump*2*pi;
    unph(:,:,:,echo) = unph(:,:,:,echo);
end

nii = make_nii(unph,voxel_size);
save_nii(nii,'unph_bestpath_corr.nii');


% fit phase images with echo times
disp('--> magnitude weighted LS fit of phase to TE ...');
[tfs0, fit_residual0] = echofit(unph,mag,TE,0); 
nii = make_nii(tfs0,voxel_size);
save_nii(nii,'tfs0.nii');
nii = make_nii(fit_residual0,voxel_size);
save_nii(nii,'fit_residual0.nii');

% [tfs1, fit_residual1] = echofit(unph,mag,TE,1); 
% nii = make_nii(tfs1,voxel_size);
% save_nii(nii,'tfs1.nii');
% nii = make_nii(fit_residual1,voxel_size);
% save_nii(nii,'fit_residual1.nii');

r_mask = 1;

for fit_thr = [20, 40]
    % extra filtering according to fitting residuals
    if r_mask
        % generate reliability map
        fit_residual_blur = smooth3(fit_residual0,'box',round(0.1./voxel_size)*2+1); 
        nii = make_nii(fit_residual_blur,voxel_size);
        save_nii(nii,'fit_residual_blur0.nii');
        R = ones(size(fit_residual_blur));
        R(fit_residual_blur >= fit_thr) = 0;
    else
        R = 1;
    end


    % normalize to main field
    % ph = gamma*dB*TE
    % dB/B = ph/(gamma*TE*B0)
    % units: TE s, gamma 2.675e8 rad/(sT), B0 3T
    tfs = -tfs0/(CF)*1e6; % unit ppm

    nii = make_nii(tfs,voxel_size);
    save_nii(nii,'tfs.nii');


    % disp('--> RESHARP to remove background field ...');
    % smv_rad = 0.3;
    % tik_reg = 5e-4;
    % % tik_reg = 0;
    % cgs_num = 500;
    % tv_reg = 2e-4;
    % z_prjs = B0_dir;
    % inv_num = 500;

    % [lfs_resharp, mask_resharp] = resharp(tfs,mask.*R,voxel_size,smv_rad,tik_reg,cgs_num);
    % % 3D 2nd order polyfit to remove any residual background
    % lfs_resharp= lfs_resharp - poly3d(lfs_resharp,mask_resharp);

    % % save nifti
    % [~,~,~] = mkdir('RESHARP');
    % nii = make_nii(lfs_resharp,voxel_size);
    % save_nii(nii,['RESHARP/lfs_resharp0_tik_', num2str(tik_reg), '_num_', num2str(cgs_num), '_fit_thr' num2str(fit_thr) '.nii']);

    % % % inversion of susceptibility 
    % % disp('--> TV susceptibility inversion on RESHARP...');
    % % sus_resharp = tvdi(lfs_resharp,mask_resharp,voxel_size,tv_reg,mag(:,:,:,end),z_prjs,inv_num); 
    % % % save nifti
    % % nii = make_nii(sus_resharp.*mask_resharp,voxel_size);
    % % save_nii(nii,['RESHARP/sus_resharp_tik_', num2str(tik_reg), '_tv_', num2str(tv_reg), '_num_', num2str(inv_num), '.nii']);

    % % iLSQR method
    % chi_iLSQR = QSM_iLSQR(lfs_resharp*(CF)/1e6,mask_resharp,'H',z_prjs,'voxelsize',voxel_size,'niter',50,'TE',1000,'B0',9.4);
    % nii = make_nii(chi_iLSQR,voxel_size);
    % save_nii(nii,['RESHARP/chi_iLSQR_resharp_fit_thr' num2str(fit_thr) '.nii']);


    % % RESHARP + MEDI
    % %%%%% normalize signal intensity by noise to get SNR %%%
    % %%%% Generate the Magnitude image %%%%
    % iMag = sqrt(sum(mag.^2,4));
    % iFreq = [];
    % N_std = 1;

    % %RDF = lfs_resharp*CF*delta_TE*1e-6;
    % RDF = lfs_resharp*CF*2*pi*delta_TE*1e-6;
    % Mask = mask_resharp;
    % save RDF.mat RDF iFreq iMag N_std Mask matrix_size...
    %      voxel_size delta_TE CF B0_dir;
    % QSM = MEDI_L1('lambda',10000);
    % nii = make_nii(QSM.*Mask,voxel_size);
    % save_nii(nii,['RESHARP/MEDI10000_RESHARP_smvrad' num2str(smv_rad) '_fit_thr' num2str(fit_thr) '.nii']);



    % V-SHARP + iLSQR
    mkdir VSHARP
    B0 =9.4;
    voxelsize = voxel_size;
    padsize = [12 12 12];
    smvsize = 12;
    [TissuePhase3d, mask_vsharp] = V_SHARP(tfs ,single(mask.*R),'smvsize',smvsize,'voxelsize',voxelsize*10);
    nii = make_nii(TissuePhase3d,voxel_size);
    save_nii(nii,['VSHARP/VSHARP_fit_thr' num2str(fit_thr) '.nii']);

    chi_iLSQR_vsharp = QSM_iLSQR(TissuePhase3d*(CF)/1e6,mask_vsharp,'H',z_prjs,'voxelsize',voxel_size,'niter',50,'TE',1000,'B0',9.4);
    nii = make_nii(chi_iLSQR_vsharp,voxel_size);
    save_nii(nii,['VSHARP/chi_iLSQR_vsharp_fit_thr' num2str(fit_thr) '.nii']);
end % end fit thr loop


% % TFI method
%     disp('TFI');
%     mkdir TFI
%     %%%% Generate the Magnitude image %%%%
%     iMag = sqrt(sum(mag.^2,4));
%     matrix_size = single(imsize(1:3));
%     voxel_size = voxelsize;
%     delta_TE = TE(2) - TE(1);
%     B0_dir = z_prjs';
%     CF = 42.6036*9.4 *1e6;
%     N_std = 1;

%     % (4) TFI of 3 voxels erosion
%     iFreq = tfs*CF*delta_TE*1e-6;
%     % erode the mask (full mask to 3mm erosion)
%     % apply R
%     mask = mask.*R;
%     % mask_erosion
%     r = 1; 
%     [X,Y,Z] = ndgrid(-r:r,-r:r,-r:r);
%     h = (X.^2/r^2 + Y.^2/r^2 + Z.^2/r^2 <= 1);
%     ker = h/sum(h(:));
%     imsize = size(mask);
%     mask_tmp = convn(mask,ker,'same');
%     mask_ero1 = zeros(imsize);
%     mask_ero1(mask_tmp > 0.999999) = 1; % no error tolerance
%     r = 2; 
%     [X,Y,Z] = ndgrid(-r:r,-r:r,-r:r);
%     h = (X.^2/r^2 + Y.^2/r^2 + Z.^2/r^2 <= 1);
%     ker = h/sum(h(:));
%     imsize = size(mask);
%     mask_tmp = convn(mask,ker,'same');
%     mask_ero2 = zeros(imsize);
%     mask_ero2(mask_tmp > 0.999999) = 1; % no error tolerance
%     r = 3; 
%     [X,Y,Z] = ndgrid(-r:r,-r:r,-r:r);
%     h = (X.^2/r^2 + Y.^2/r^2 + Z.^2/r^2 <= 1);
%     ker = h/sum(h(:));
%     imsize = size(mask);
%     mask_tmp = convn(mask,ker,'same');
%     mask_ero3 = zeros(imsize);
%     mask_ero3(mask_tmp > 0.999999) = 1; % no error tolerance


%     Mask = mask;
%     Mask_G = Mask;
%     P_B = 30;
%     P = 1 * Mask + P_B * (1-Mask);
%     RDF = 0;
%     save RDF_brain.mat matrix_size voxel_size delta_TE B0_dir CF iMag N_std iFreq Mask Mask_G P RDF
%     % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 500);
%     % nii = make_nii(QSM.*Mask,vox);
%     % save_nii(nii,'TFI/TFI_lambda500_full.nii');
%     QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 20000);
%     nii = make_nii(QSM.*Mask,voxelsize);
%     save_nii(nii,'TFI/TFI_lambda20000_full.nii');
%     % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1500);
%     % nii = make_nii(QSM.*Mask,vox);
%     % save_nii(nii,'TFI/TFI_lambda1500_full.nii');
%     % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 2000);
%     % nii = make_nii(QSM.*Mask,vox);
%     % save_nii(nii,'TFI/TFI_lambda2000_full.nii');


%     Mask = mask_ero1;
%     Mask_G = Mask;
%     P_B = 30;
%     P = 1 * Mask + P_B * (1-Mask);
%     RDF = 0;
%     save RDF_brain.mat matrix_size voxel_size delta_TE B0_dir CF iMag N_std iFreq Mask Mask_G P RDF
%     % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 500);
%     % nii = make_nii(QSM.*Mask,vox);
%     % save_nii(nii,'TFI/TFI_lambda500_ero1.nii');
%     QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 20000);
%     nii = make_nii(QSM.*Mask,voxelsize);
%     save_nii(nii,'TFI/TFI_lambda20000_ero1.nii');
%     % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1500);
%     % nii = make_nii(QSM.*Mask,vox);
%     % save_nii(nii,'TFI/TFI_lambda1500_ero1.nii');
%     % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 2000);
%     % nii = make_nii(QSM.*Mask,vox);
%     % save_nii(nii,'TFI/TFI_lambda2000_ero1.nii');


%     Mask = mask_ero2;
%     Mask_G = Mask;
%     P_B = 30;
%     P = 1 * Mask + P_B * (1-Mask);
%     RDF = 0;
%     save RDF_brain.mat matrix_size voxel_size delta_TE B0_dir CF iMag N_std iFreq Mask Mask_G P RDF
%     % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 500);
%     % nii = make_nii(QSM.*Mask,vox);
%     % save_nii(nii,'TFI/TFI_lambda500_ero2.nii');
%     QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 20000);
%     nii = make_nii(QSM.*Mask,voxelsize);
%     save_nii(nii,'TFI/TFI_lambda20000_ero2.nii');
%     % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1500);
%     % nii = make_nii(QSM.*Mask,vox);
%     % save_nii(nii,'TFI/TFI_lambda1500_ero2.nii');
%     % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 2000);
%     % nii = make_nii(QSM.*Mask,vox);
%     % save_nii(nii,'TFI/TFI_lambda2000_ero2.nii');


%     Mask = mask_ero3;
%     Mask_G = Mask;
%     P_B = 30;
%     P = 1 * Mask + P_B * (1-Mask);
%     RDF = 0;
%     save RDF_brain.mat matrix_size voxel_size delta_TE B0_dir CF iMag N_std iFreq Mask Mask_G P RDF
%     % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 500);
%     % nii = make_nii(QSM.*Mask,vox);
%     % save_nii(nii,'TFI/TFI_lambda500_ero3.nii');
%     QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 20000);
%     nii = make_nii(QSM.*Mask,voxelsize);
%     save_nii(nii,'TFI/TFI_lambda20000_ero3.nii');
%     % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1500);
%     % nii = make_nii(QSM.*Mask,vox);
%     % save_nii(nii,'TFI/TFI_lambda1500_ero3.nii');
%     % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 2000);
%     % nii = make_nii(QSM.*Mask,vox);
%     % save_nii(nii,'TFI/TFI_lambda2000_ero3.nii');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[iFreq_raw N_std] = Fit_ppm_complex(iField);
nii = make_nii(iFreq_raw,voxel_size);
save_nii(nii,'iFreq_raw.nii');

% phase unwrapping using prelude
!prelude -a src/corr_mag1.nii -p iFreq_raw.nii -u iFreq_un -m mask.nii
!gunzip iFreq_un.nii.gz
nii = load_nii('iFreq_un.nii');
iFreq = double(nii.img);


% % RESHARP background field removal
% [RDF, mask_resharp] = resharp(iFreq,mask,voxel_size,0.3,5e-4,500);
% nii = make_nii(RDF,voxel_size);
% save_nii(nii,'RDF_resharp.nii');

% lfs_resharp = RDF/(CF*delta_TE*1e-6);


% % iLSQR method
% chi_iLSQR_0 = QSM_iLSQR(lfs_resharp*(CF)/1e6,mask_resharp,'H',z_prjs,'voxelsize',voxel_size,'niter',50,'TE',1000,'B0',9.4);
% nii = make_nii(chi_iLSQR_0,voxel_size);
% save_nii(nii,'chi_iLSQR_0.nii');

tfs = iFreq/(CF*delta_TE*1e-6);

% V-SHARP + iLSQR
B0 =9.4;
voxelsize = voxel_size;
padsize = [12 12 12];
smvsize = 12;
[TissuePhase3d, mask_vsharp] = V_SHARP(tfs ,single(mask),'smvsize',smvsize,'voxelsize',voxelsize*10);
nii = make_nii(TissuePhase3d,voxel_size);
save_nii(nii,'VSHARP.nii');

chi_iLSQR_0 = QSM_iLSQR(TissuePhase3d*(CF)/1e6,mask_vsharp,'H',z_prjs,'voxelsize',voxel_size,'niter',50,'TE',1000,'B0',9.4);
nii = make_nii(chi_iLSQR_0,voxel_size);
save_nii(nii,'chi_cpx_iLSQR_0_vsharp.nii');


% % TVDI
% sus_resharp = tvdi(lfs_resharp,mask_resharp,voxel_size,5e-5,abs(iField(:,:,:,end)),B0_dir,500);
% nii = make_nii(sus_resharp.*mask_resharp,voxel_size);
% save_nii(nii,'QSM_resharp_5e-5.nii');

% sus_resharp = tvdi(lfs_resharp,mask_resharp,voxel_size,1e-5,abs(iField(:,:,:,end)),B0_dir,500);
% nii = make_nii(sus_resharp.*mask_resharp,voxel_size);
% save_nii(nii,'QSM_resharp_1e-5.nii');

% sus_resharp = tvdi(lfs_resharp,mask_resharp,voxel_size,2e-4,abs(iField(:,:,:,end)),B0_dir,500);
% nii = make_nii(sus_resharp.*mask_resharp,voxel_size);
% save_nii(nii,'QSM_resharp_2e-4.nii');

% sus_resharp = tvdi(lfs_resharp,mask_resharp,voxel_size,1e-3,abs(iField(:,:,:,end)),B0_dir,500);
% nii = make_nii(sus_resharp.*mask_resharp,voxel_size);
% save_nii(nii,'QSM_resharp_1e-3.nii');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mkdir TFS_TIK_PRE_ERO0
% cd TFS_TIK_PRE_ERO0
% tfs_pad = padarray(iFreq/(CF*delta_TE*1e-6),[0 0 20]);
% mask_pad = padarray(mask,[0 0 20]);
% % R_pad = padarray(R,[0 0 20]);
% r=0;
% Tik_weight = [1e-3];
% TV_weight = 2e-4;
% for i = 1:length(Tik_weight)
% 	% chi = tikhonov_qsm(tfs_pad, mask_pad, 1, mask_pad, mask_pad, TV_weight, Tik_weight(i), vox, z_prjs, 200);
% 	% nii = make_nii(chi(:,:,21:end-20).*mask_pad(:,:,21:end-20),vox);
% 	% save_nii(nii,['TIK_ero' num2str(r) '_TV_' num2str(TV_weight) '_Tik_' num2str(Tik_weight(i)) '_PRE_200.nii']);
% 	% chi = tikhonov_qsm(tfs_pad, mask_pad, 1, mask_pad, mask_pad, TV_weight, Tik_weight(i), vox, z_prjs, 500);
% 	% nii = make_nii(chi(:,:,21:end-20).*mask_pad(:,:,21:end-20),vox);
% 	% save_nii(nii,['TIK_ero' num2str(r) '_TV_' num2str(TV_weight) '_Tik_' num2str(Tik_weight(i)) '_PRE_500.nii']);
% 	% chi = tikhonov_qsm(tfs_pad, mask_pad, 1, mask_pad, mask_pad, TV_weight, Tik_weight(i), vox, z_prjs, 2000);
% 	% nii = make_nii(chi(:,:,21:end-20).*mask_pad(:,:,21:end-20),vox);
% 	% save_nii(nii,['TIK_ero' num2str(r) '_TV_' num2str(TV_weight) '_Tik_' num2str(Tik_weight(i)) '_PRE_2000.nii']);
% 	chi = tikhonov_qsm(tfs_pad, mask_pad, 1, mask_pad, mask_pad, TV_weight, Tik_weight(i), voxel_size, B0_dir, 2000);
% 	nii = make_nii(chi(:,:,21:end-20).*mask_pad(:,:,21:end-20),voxel_size);
% 	save_nii(nii,['TIK_ero' num2str(r) '_TV_' num2str(TV_weight) '_Tik_' num2str(Tik_weight(i)) '_PRE_2000.nii']);
% end
% cd ..


cd ..
