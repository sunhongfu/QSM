% read in combined magnitude and phase images
path_mag = '/home/hongfu/NCIgb5_scratch/Jon_temp_downloads/point75_MP2RAGE_Multiechos/04_BH_451/MAG_UNCOMB';
path_ph = '/home/hongfu/NCIgb5_scratch/Jon_temp_downloads/point75_MP2RAGE_Multiechos/04_BH_451/PHA_UNCOMB';
path_out = '/home/hongfu/NCIgb5_scratch/hongfu/ME-MP2RAGE/04_BH_451';

%% read in DICOMs of both uncombined magnitude and raw unfiltered phase images
path_mag = cd(cd(path_mag));
mag_list = dir([path_mag '/*.dcm']);
mag_list = mag_list(~strncmpi('.', {mag_list.name}, 1));
path_ph = cd(cd(path_ph));
ph_list = dir([path_ph '/*.dcm']);
ph_list = ph_list(~strncmpi('.', {ph_list.name}, 1));

% number of total (2 inversion) slices (mag and ph should be the same)
nSL = length(ph_list);


% get the sequence parameters
dicom_info = dicominfo([path_ph,filesep,ph_list(end).name]);
NumberOfEchoes = dicom_info.EchoNumber; 

% get the sequence parameters
for i = 1:nSL/2/NumberOfEchoes:nSL/2 % read in TEs
    dicom_info = dicominfo([path_ph,filesep,ph_list(i).name]);
    TE(dicom_info.EchoNumber) = dicom_info.EchoTime*1e-3;
end
vox = [dicom_info.PixelSpacing(1), dicom_info.PixelSpacing(2), dicom_info.SliceThickness];

% angles (z projections of the image x y z coordinates) 
Xz = dicom_info.ImageOrientationPatient(3);
Yz = dicom_info.ImageOrientationPatient(6);
Zxyz = cross(dicom_info.ImageOrientationPatient(1:3),dicom_info.ImageOrientationPatient(4:6));
Zz = Zxyz(3);
z_prjs = [Xz, Yz, Zz];

% read in measurements
mag = zeros(dicom_info.Rows,dicom_info.Columns,nSL,'single');
ph = zeros(dicom_info.Rows,dicom_info.Columns,nSL,'single');
for i = 1:nSL
    mag(:,:,i) = single(dicomread([path_mag,filesep,mag_list(i).name]));
    ph(:,:,i) = single(dicomread([path_ph,filesep,ph_list(i).name]));
end


% crop mosaic into individual images
AcqMatrix = regexp(dicom_info.Private_0051_100b,'(\d)*(\d)','match');
if strcmpi(dicom_info.InPlanePhaseEncodingDirection,'COL') % A/P
% phase encoding along column
    wRow = round(str2num(AcqMatrix{1})/dicom_info.PercentSampling*100);
    wCol = str2num(AcqMatrix{2});
else % L/R
    wCol = round(str2num(AcqMatrix{1})/dicom_info.PercentSampling*100);
    wRow = str2num(AcqMatrix{2});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% works for square 32 channels (faster)
mag = mat2cell(mag,[wRow wRow wRow wRow wRow wRow], [wCol wCol wCol wCol wCol wCol], nSL);
mag_all = cat(4,mag{1,1}, mag{1,2}, mag{1,3}, mag{1,4}, mag{1,5}, mag{1,6}, mag{2,1}, mag{2,2}, mag{2,3}, mag{2,4}, mag{2,5}, mag{2,6}, mag{3,1}, mag{3,2}, mag{3,3}, mag{3,4}, mag{3,5}, mag{3,6}, mag{4,1}, mag{4,2}, mag{4,3}, mag{4,4}, mag{4,5}, mag{4,6}, mag{5,1}, mag{5,2}, mag{5,3}, mag{5,4}, mag{5,5}, mag{5,6}, mag{6,1}, mag{6,2});
clear mag
mag_all = reshape(mag_all, wRow, wCol, nSL/2/NumberOfEchoes, NumberOfEchoes, 2, 32);
mag_all = permute(mag_all,[2 1 3 4 6 5]);

ph = mat2cell(ph,[wRow wRow wRow wRow wRow wRow], [wCol wCol wCol wCol wCol wCol], nSL);
ph_all = cat(4,ph{1,1}, ph{1,2}, ph{1,3}, ph{1,4}, ph{1,5}, ph{1,6}, ph{2,1}, ph{2,2}, ph{2,3}, ph{2,4}, ph{2,5}, ph{2,6}, ph{3,1}, ph{3,2}, ph{3,3}, ph{3,4}, ph{3,5}, ph{3,6}, ph{4,1}, ph{4,2}, ph{4,3}, ph{4,4}, ph{4,5}, ph{4,6}, ph{5,1}, ph{5,2}, ph{5,3}, ph{5,4}, ph{5,5}, ph{5,6}, ph{6,1}, ph{6,2});
clear ph
ph_all = reshape(ph_all, wRow, wCol, nSL/2/NumberOfEchoes, NumberOfEchoes, 2, 32);
ph_all = permute(ph_all,[2 1 3 4 6 5]);
ph_all = 2*pi.*(ph_all - single(dicom_info.SmallestImagePixelValue))/(single(dicom_info.LargestImagePixelValue - dicom_info.SmallestImagePixelValue)) - pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % uncombined mosiac
% % read in measurements
% mag = zeros(dicom_info.Rows,dicom_info.Columns,nSL,'single');
% ph = zeros(dicom_info.Rows,dicom_info.Columns,nSL,'single');
% for i = 1:nSL
%     mag(:,:,i) = single(dicomread([path_mag,filesep,mag_list(i).name]));
%     ph(:,:,i) = single(dicomread([path_ph,filesep,ph_list(i).name]));
% end

% % crop mosaic into individual images
% AcqMatrix = regexp(dicom_info.Private_0051_100b,'(\d)*(\d)','match');
% if strcmpi(dicom_info.InPlanePhaseEncodingDirection,'COL') % A/P
% % phase encoding along column
%     wRow = round(str2num(AcqMatrix{1})/dicom_info.PercentSampling*100);
%     wCol = str2num(AcqMatrix{2});
% else % L/R
%     wCol = round(str2num(AcqMatrix{1})/dicom_info.PercentSampling*100);
%     wRow = str2num(AcqMatrix{2});
% end

% nCol = double(dicom_info.Columns/wCol);
% nRow = double(dicom_info.Rows/wRow);
% nChan = double(dicom_info.Private_0019_100a);

% mag_all = zeros(wRow,wCol,nChan,nSL,'single');
% ph_all = zeros(wRow,wCol,nChan,nSL,'single');
% for i = 1:nSL
%     for x = 1:wRow
%         for y = 1:wCol
%             for z = 1:nChan
%                 X = floor((z-1)/nCol)*wRow + x;
%                 Y = mod(z-1,nCol)*wCol + y;
%                 mag_all(x,y,z,i) = mag(X,Y,i);
%                 ph_all(x,y,z,i) = ph(X,Y,i);
%             end
%         end
%     end
% end


% % reshape and permute into COLS, ROWS, SLICES, ECHOES, CHANS
% mag_all = reshape(mag_all,[wRow,wCol,nChan,nSL/2/NumberOfEchoes,NumberOfEchoes,2]);
% mag_all = permute(mag_all,[2 1 4 5 3 6]);
% ph_all = reshape(ph_all,[wRow,wCol,nChan,nSL/2/NumberOfEchoes,NumberOfEchoes,2]);
% ph_all = permute(ph_all,[2 1 4 5 3 6]);
% % 0028,0106  Smallest Image Pixel Value: 0
% % 0028,0107  Largest Image Pixel Value: 4094
% % conver scale to -pi to pi
% ph_all = 2*pi.*(ph_all - single(dicom_info.SmallestImagePixelValue))./(single(dicom_info.LargestImagePixelValue - dicom_info.SmallestImagePixelValue)) - pi;

% imsize = size(ph_all);

% clear mag ph

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define output directories
path_qsm = [path_out '/QSM_MEMP2RAGE_7T'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use only the second inversion
mag_all_i2 = mag_all(:,:,:,:,:,2);
ph_all_i2 = ph_all(:,:,:,:,:,2);
clear mag_all ph_all
imsize = size(ph_all_i2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save the raw data for future use
save('raw.mat','-v7.3');


% BEGIN THE QSM RECON PIPELINE
% initial quick brain mask
% simple sum-of-square combination
mag1_sos = sqrt(sum(mag_all_i2(:,:,:,1,:).^2,5));
nii = make_nii(mag1_sos,vox);
save_nii(nii,'mag1_sos.nii');

unix('N4BiasFieldCorrection -i mag1_sos.nii -o mag1_sos_n4.nii');

unix('bet2 mag1_sos_n4.nii BET -f 0.3 -m');
% set a lower threshold for postmortem
% unix('bet2 mag1_sos.nii BET -f 0.1 -m');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);

% coil combination % smoothing factor 10?
ph_corr = zeros(imsize(1:4));
mag_corr = zeros(imsize(1:4));
% [ph_corr,mag_corr] = geme_cmb(mag_all_i2.*exp(1j*ph_all_i2),vox,TE,mask);

[ph_corr(:,:,:,1:2:end),mag_corr(:,:,:,1:2:end)] = geme_cmb(mag_all_i2(:,:,:,1:2:end,:).*exp(1j*ph_all_i2(:,:,:,1:2:end,:)),vox,TE(1:2:end),mask,'gaussian');
!mkdir odd_gaussian; mv offsets* odd_gaussian
[ph_corr(:,:,:,2:2:end),mag_corr(:,:,:,2:2:end)] = geme_cmb(mag_all_i2(:,:,:,2:2:end,:).*exp(1j*ph_all_i2(:,:,:,2:2:end,:)),vox,TE(2:2:end),mask,'gaussian');
!mkdir even_gaussian; mv offsets* even_gaussian

clear mag_all_i2 ph_all_i2

% save niftis after coil combination
mkdir('src');
for echo = 1:imsize(4)
    nii = make_nii(mag_corr(:,:,:,echo),vox);
    save_nii(nii,['src/mag_corr' num2str(echo) '.nii']);
    nii = make_nii(ph_corr(:,:,:,echo),vox);
    save_nii(nii,['src/ph_corr' num2str(echo) '.nii']);
end

save('raw.mat','ph_corr','mag_corr','mask','-append');




% (1) unwrap the phase using best path
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

unph = zeros(imsize(1:4));

for echo_num = 1:imsize(4)
    setenv('echo_num',num2str(echo_num));
    fid = fopen(['wrapped_phase' num2str(echo_num) '.dat'],'w');
    fwrite(fid,ph_corr(:,:,:,echo_num),'float');
    fclose(fid);

    bash_script = ['${pathstr}/3DSRNCP wrapped_phase${echo_num}.dat mask_unwrp.dat ' ...
        'unwrapped_phase${echo_num}.dat $nv $np $ns reliability${echo_num}.dat'];
    unix(bash_script) ;

    fid = fopen(['unwrapped_phase' num2str(echo_num) '.dat'],'r');
    tmp = fread(fid,'float');
    % tmp = tmp - tmp(1);
    unph(:,:,:,echo_num) = reshape(tmp - round(mean(tmp(mask==1))/(2*pi))*2*pi ,imsize(1:3)).*mask;
    fclose(fid);
end

nii = make_nii(unph,vox);
save_nii(nii,'unph_bestpath_before_jump_correction.nii');

% remove all the temp files
! rm *.dat

% 2pi jumps correction
nii = load_nii('unph_diff.nii');
unph_diff = double(nii.img);
unph_diff = unph_diff/2;
for echo = 2:imsize(4)
    meandiff = unph(:,:,:,echo)-unph(:,:,:,1)-double(echo-1)*unph_diff;
    meandiff = meandiff(mask==1);
    meandiff = mean(meandiff(:));
    njump = round(meandiff/(2*pi));
    disp(['    ' num2str(njump) ' 2pi jumps for TE' num2str(echo)]);
    unph(:,:,:,echo) = unph(:,:,:,echo) - njump*2*pi;
    unph(:,:,:,echo) = unph(:,:,:,echo).*mask;
end
nii = make_nii(unph,vox);
save_nii(nii,'unph_bestpath.nii');

save('raw.mat','unph','-append');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set parameters
fit_thr = 5;
tik_reg = 1e-6;
cgs_num = 500;
lsqr_num = 500;
% tv_reg = 2e-4;
inv_num = 500;


% fit phase images with echo times
disp('--> magnitude weighted LS fit of phase to TE ...');
[tfs_0, fit_residual_0] = echofit(unph,mag_corr,TE,0); 
% normalize to main field
% ph = gamma*dB*TE
% dB/B = ph/(gamma*TE*B0)
% units: TE s, gamma 2.675e8 rad/(sT), B0 7T
tfs_0 = tfs_0/(2.675e8*dicom_info.MagneticFieldStrength)*1e6; % unit ppm
nii = make_nii(tfs_0,vox);
save_nii(nii,'tfs_0.nii');
nii = make_nii(fit_residual_0,vox);
save_nii(nii,'fit_residual_0.nii');
% extra filtering according to fitting residuals
% generate reliability map
fit_residual_0_blur = smooth3(fit_residual_0,'box',round(1./vox)*2+1); 
nii = make_nii(fit_residual_0_blur,vox);
save_nii(nii,'fit_residual_0_blur.nii');
R_0 = ones(size(fit_residual_0_blur));
R_0(fit_residual_0_blur >= fit_thr) = 0;

save('raw.mat','tfs_0','fit_residual_0','fit_residual_0_blur','R_0','fit_thr','cgs_num','tik_reg','inv_num','-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('raw.mat','tfs_0','mask','R_0','vox','cgs_num','dicom_info','z_prjs','mag_corr','imsize','TE');

for smv_rad = [2]
% for smv_rad = [1 2 3]
    % RE-SHARP (tik_reg: Tikhonov regularization parameter)
    disp('--> RESHARP to remove background field ...');
    % [lfs_resharp_0, mask_resharp_0] = resharp_lsqr(tfs_0,mask.*R_0,vox,smv_rad,lsqr_num);
    [lfs_resharp_0, mask_resharp_0] = resharp(tfs_0,mask.*R_0,vox,smv_rad,tik_reg,cgs_num);
    % save nifti
    mkdir('RESHARP');
    nii = make_nii(lfs_resharp_0,vox);
    % save_nii(nii,['RESHARP/lfs_resharp_0_smvrad' num2str(smv_rad) '_lsqr.nii']);
    save_nii(nii,['RESHARP/lfs_resharp_0_smvrad' num2str(smv_rad) '_cgs_' num2str(tik_reg) '.nii']);

    % iLSQR
    chi_iLSQR_0 = QSM_iLSQR(lfs_resharp_0*(2.675e8*dicom_info.MagneticFieldStrength)/1e6,mask_resharp_0,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
    nii = make_nii(chi_iLSQR_0,vox);
    save_nii(nii,['RESHARP/chi_iLSQR_smvrad' num2str(smv_rad) '.nii']);

    % MEDI
    %%%%% normalize signal intensity by noise to get SNR %%%
    %%%% Generate the Magnitude image %%%%
    iMag = sqrt(sum(mag_corr.^2,4));
    % [iFreq_raw N_std] = Fit_ppm_complex(ph_corr);
    matrix_size = single(imsize(1:3));
    voxel_size = vox;
    delta_TE = TE(2) - TE(1);
    B0_dir = z_prjs';
    CF = dicom_info.ImagingFrequency *1e6;
    iFreq = [];
    N_std = 1;
    RDF = lfs_resharp_0*2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6;
    Mask = mask_resharp_0;
    save RDF.mat RDF iFreq iMag N_std Mask matrix_size...
         voxel_size delta_TE CF B0_dir;
    QSM = MEDI_L1('lambda',1000);
    nii = make_nii(QSM.*Mask,vox);
    save_nii(nii,['RESHARP/MEDI1000_RESHARP_smvrad' num2str(smv_rad) '.nii']);
    
    % %TVDI
    % sus_resharp = tvdi(lfs_resharp_0,mask_resharp_0,vox,2e-4,iMag,z_prjs,500); 
    % nii = make_nii(sus_resharp.*mask_resharp_0,vox);
    % save_nii(nii,['RESHARP/TV_2e-4_smvrad' num2str(smv_rad) '.nii']);

end


mkdir('VSHARP');
% V-SHARP + iLSQR
padsize = [12 12 12];
smvsize = 12;
[lfs_vsharp, mask_vsharp] = V_SHARP(tfs_0 ,single(mask.*R_0),'smvsize',smvsize,'voxelsize',vox);
nii = make_nii(lfs_vsharp,vox);
save_nii(nii,'VSHARP/lfs_vsharp.nii');

chi_iLSQR_0 = QSM_iLSQR(lfs_vsharp*(2.675e8*dicom_info.MagneticFieldStrength)/1e6,mask_vsharp,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
nii = make_nii(chi_iLSQR_0,vox);
save_nii(nii,'VSHARP/chi_iLSQR_0_vsharp_iter50.nii');

% MEDI
%%%%% normalize signal intensity by noise to get SNR %%%
%%%% Generate the Magnitude image %%%%
iMag = sqrt(sum(mag_corr.^2,4));
% [iFreq_raw N_std] = Fit_ppm_complex(ph_corr);
matrix_size = single(imsize(1:3));
voxel_size = vox;
delta_TE = TE(2) - TE(1);
B0_dir = z_prjs';
CF = dicom_info.ImagingFrequency *1e6;
iFreq = [];
N_std = 1;
RDF = lfs_vsharp*2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6;
Mask = mask_vsharp;
save RDF.mat RDF iFreq iMag N_std Mask matrix_size...
     voxel_size delta_TE CF B0_dir;
QSM = MEDI_L1('lambda',1000);
nii = make_nii(QSM.*Mask,vox);
save_nii(nii,['VSHARP/MEDI1000_RESHARP_smvrad' num2str(smv_rad) '.nii']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% perform inversion on each individual echo
tfs_resharp = zeros(imsize(1:4));
chi_iLSQR = zeros(imsize(1:4));
smv_rad = 2;
for i = 1:imsize(4)
    % RE-SHARP (tik_reg: Tikhonov regularization parameter)
    disp('--> RESHARP to remove background field ...');
    [tfs_resharp(:,:,:,i), mask_resharp] = resharp(unph(:,:,:,i),mask.*R_0,vox,smv_rad,tik_reg,cgs_num);
    nii = make_nii(tfs_resharp(:,:,:,i),vox);
    save_nii(nii,['RESHARP/RESHARP_' 'e' num2str(i) '.nii']);
    % iLSQR
    chi_iLSQR(:,:,:,i) = QSM_iLSQR(tfs_resharp(:,:,:,i),mask_resharp,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000*TE(i),'B0',dicom_info.MagneticFieldStrength);
    nii = make_nii(chi_iLSQR(:,:,:,i),vox);
    save_nii(nii,['RESHARP/chi_iLSQR_50_' 'e' num2str(i) '.nii']);
end
chi_iLSQR_ave = mean(chi_iLSQR,4);
nii = make_nii(chi_iLSQR_ave,vox);
save_nii(nii,'RESHARP/chi_iLSQR_ave.nii');



% perform inversion on each individual echo
tfs_vsharp = zeros(imsize(1:4));
chi_iLSQR = zeros(imsize(1:4));
padsize = [12 12 12];
smvsize = 12;
for i = 1:imsize(4)
    % VSHARP (tik_reg: Tikhonov regularization parameter)
    disp('--> VSHARP to remove background field ...');
    [tfs_vsharp(:,:,:,i), mask_vsharp] = V_SHARP(unph(:,:,:,i) ,single(mask.*R_0),'smvsize',smvsize,'voxelsize',vox);
    nii = make_nii(tfs_vsharp(:,:,:,i),vox);
    save_nii(nii,['VSHARP/VSHARP_' 'e' num2str(i) '.nii']);
    % iLSQR
    chi_iLSQR(:,:,:,i) = QSM_iLSQR(tfs_vsharp(:,:,:,i),mask_vsharp,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000*TE(i),'B0',dicom_info.MagneticFieldStrength);
    nii = make_nii(chi_iLSQR(:,:,:,i),vox);
    save_nii(nii,['VSHARP/chi_iLSQR_50_' 'e' num2str(i) '.nii']);
end
chi_iLSQR_ave = mean(chi_iLSQR,4);
nii = make_nii(chi_iLSQR_ave,vox);
save_nii(nii,'VSHARP/chi_iLSQR_ave.nii');







% %% only using the last echo
% ph_last = ph(:,:,:,end);
% mag_last = mag(:,:,:,end);

% laplacian unwrapping of the last echo
disp('--> unwrap aliasing phase using laplacian...');
Options.voxelSize = vox;
unph_lap = lapunwrap(iFreq_raw, Options);
nii = make_nii(unph_lap, vox);
save_nii(nii,'unph_lap.nii');

% unph_lap = LaplacianPhaseUnwrap(iFreq_raw,'voxelsize',vox);
[TissuePhase3d, mask_vsharp] = V_SHARP(-unph_lap ,single(mask),'smvsize',12,'voxelsize',vox*10);
nii = make_nii(TissuePhase3d,vox);
save_nii(nii,['VSHARP.nii']);

    chi_iLSQR_0 = QSM_iLSQR(TissuePhase3d,mask_vsharp,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',9.4);
    nii = make_nii(chi_iLSQR_0,voxel_size);
    save_nii(nii,['chi_iLSQR_0_vsharp' num2str(echo) '.nii']);

 


% complex fitting
iField = mag.*exp(1i*ph);
%bipolar correction first
iField = iField_correction(iField,vox);
[iFreq_raw N_std] = Fit_ppm_complex(iField);
nii = make_nii(iFreq_raw,vox);
save_nii(nii,'iFreq_raw.nii');
nii = make_nii(N_std,vox);
save_nii(nii,'N_std.nii');


% complex fitting
iField_corr = mag_corr.*exp(1i*ph_corr);
iField_corr = iField_correction(iField_corr,vox);
[iFreq_raw_corr N_std_corr] = Fit_ppm_complex(iField_corr);
nii = make_nii(iFreq_raw_corr,vox);
save_nii(nii,'iFreq_raw_corr.nii');
nii = make_nii(N_std_corr,vox);
save_nii(nii,'N_std_corr.nii');


% laplacian unwrapping
Unwrapped_Phase = LaplacianPhaseUnwrap(iFreq_raw_corr,'voxelsize',vox);
nii = make_nii(Unwrapped_Phase,vox);
save_nii(nii,'Unwrapped_Phase.nii');


% FUDGE unwrapping
unph_fudge = fudge(iFreq_raw);
nii = make_nii(unph_fudge,vox);
save_nii(nii,'unph_fudge.nii');


unph_fudge_corr = fudge(iFreq_raw_corr);
nii = make_nii(unph_fudge_corr,vox);
save_nii(nii,'unph_fudge_corr.nii');


% V-SHARP + iLSQR
voxelsize = vox;
padsize = [12 12 12];
smvsize = 12;
[TissuePhase3d, mask_vsharp] = V_SHARP(tfs_0 ,single(mask),'smvsize',smvsize,'voxelsize',vox);
nii = make_nii(TissuePhase3d,vox);
save_nii(nii,'VSHARP.nii');

delta_TE = TE(2) - TE(1);
% tfs = TissuePhase3d/(2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6);

chi_iLSQR_0 = QSM_iLSQR(-TissuePhase3d,mask_vsharp,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',delta_TE*1000,'B0',dicom_info.MagneticFieldStrength);
nii = make_nii(chi_iLSQR_0,vox);
save_nii(nii,'chi_iLSQR_0_vsharp.nii');



%%%%%%%%%%%%%%%%%%%%%%%
% fudge + v-sharp + iLSQR for each echo
TissuePhase3d = zeros(size(mag));
mask_vsharp = zeros(size(mag));
chi_iLSQR_0 = zeros(size(mag));
unph_fudge = zeros(size(mag));

for echo = 1:4
    unph_fudge(:,:,:,echo) = fudge(ph(:,:,:,echo));
    [TissuePhase3d(:,:,:,echo), mask_vsharp(:,:,:,echo)] = V_SHARP(-unph_fudge(:,:,:,echo) ,single(mask),'smvsize',smvsize,'voxelsize',vox);
    nii = make_nii(TissuePhase3d(:,:,:,echo),vox);
    save_nii(nii,['VSHARP' num2str(echo) '.nii']);

    chi_iLSQR_0(:,:,:,echo) = QSM_iLSQR(TissuePhase3d(:,:,:,echo),mask_vsharp(:,:,:,echo),'H',z_prjs,'voxelsize',vox,'niter',50,'TE',TE(echo)*1000,'B0',dicom_info.MagneticFieldStrength);
    nii = make_nii(chi_iLSQR_0(:,:,:,echo),vox);
    save_nii(nii,['chi_iLSQR_0_vsharp' num2str(echo) '.nii']);
end

chi_iLSQR_0_ave = mean(chi_iLSQR_0,4);
nii = make_nii(chi_iLSQR_0_ave,vox);
save_nii(nii,'chi_iLSQR_0_ave.nii');


