% read in uncombined magnitude and phase images
path_mag = '/gpfs/M2Scratch/NCIgb5/hongfu/SPECIMAN/1.10.1/1.10.1.344/1.10.1.344.1.2/1.10.1.344.1.2.7/dicom_series';
path_ph = '/gpfs/M2Scratch/NCIgb5/hongfu/SPECIMAN/1.10.1/1.10.1.344/1.10.1.344.1.2/1.10.1.344.1.2.8/dicom_series';
path_out = '/home/hongfu/NCIgb5_scratch/hongfu/SPECIMAN/4';

%% read in DICOMs of both uncombined magnitude and raw unfiltered phase images

%% read in DICOMs of both uncombined magnitude and raw unfiltered phase images
path_mag = cd(cd(path_mag));
mag_list = dir([path_mag '/*.dcm']);
mag_list = mag_list(~strncmpi('.', {mag_list.name}, 1));
path_ph = cd(cd(path_ph));
ph_list = dir([path_ph '/*.dcm']);
ph_list = ph_list(~strncmpi('.', {ph_list.name}, 1));

% number of slices (mag and ph should be the same)
nSL = length(ph_list);

% get the sequence parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EchoTrainLength = 3; % this number is wrong in DICOM header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nSL/EchoTrainLength:nSL % read in TEs
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

nCol = double(dicom_info.Columns/wCol);
nRow = double(dicom_info.Rows/wRow);
nChan = double(dicom_info.Private_0019_100a);

mag_all = zeros(wRow,wCol,nChan,nSL,'single');
ph_all = zeros(wRow,wCol,nChan,nSL,'single');
for i = 1:nSL
    for x = 1:wRow
        for y = 1:wCol
            for z = 1:nChan
                X = floor((z-1)/nCol)*wRow + x;
                Y = mod(z-1,nCol)*wCol + y;
                mag_all(x,y,z,i) = mag(X,Y,i);
                ph_all(x,y,z,i) = ph(X,Y,i);
            end
        end
    end
end

% reshape and permute into COLS, ROWS, SLICES, ECHOES, CHANS
mag_all = reshape(mag_all,[wRow,wCol,nChan,nSL/EchoTrainLength,EchoTrainLength]);
mag_all = permute(mag_all,[2 1 4 5 3]);
ph_all = reshape(ph_all,[wRow,wCol,nChan,nSL/EchoTrainLength,EchoTrainLength]);
ph_all = permute(ph_all,[2 1 4 5 3]);
% 0028,0106  Smallest Image Pixel Value: 0
% 0028,0107  Largest Image Pixel Value: 4094
% conver scale to -pi to pi
ph_all = 2*pi.*(ph_all - single(dicom_info.SmallestImagePixelValue))./(single(dicom_info.LargestImagePixelValue - dicom_info.SmallestImagePixelValue)) - pi;

imsize = size(ph_all);

% define output directories
path_qsm = [path_out '/QSM_MEGE_7T'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);
% save the raw data for future use
clear mag ph
save('raw.mat','-v7.3');




% BEGIN THE QSM RECON PIPELINE
% initial quick brain mask
% simple sum-of-square combination
mag1_sos = sqrt(sum(mag_all(:,:,:,1,:).^2,5));
nii = make_nii(mag1_sos,vox);
save_nii(nii,'mag1_sos.nii');
% unix('bet2 mag1_sos.nii BET -f 0.4 -w 2 -m');
% set a lower threshold for postmortem
unix('bet2 mag1_sos.nii BET -f 0.1 -m');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);

% coil combination % smoothing factor 10?
ph_corr = zeros(imsize(1:4));
mag_corr = zeros(imsize(1:4));
[ph_corr,mag_corr] = geme_cmb(mag_all.*exp(1j*ph_all),vox,TE,mask);

% [ph_corr(:,:,:,1:2:end),mag_corr(:,:,:,1:2:end)] = geme_cmb(mag_all(:,:,:,1:2:end,:).*exp(1j*ph_all(:,:,:,1:2:end,:)),vox,TE(1:2:end),mask);
% [ph_corr(:,:,:,2:2:end),mag_corr(:,:,:,2:2:end)] = geme_cmb(mag_all(:,:,:,2:2:end,:).*exp(1j*ph_all(:,:,:,2:2:end,:)),vox,TE(2:2:end),mask);

% save niftis after coil combination
mkdir('src');
for echo = 1:imsize(4)
    nii = make_nii(mag_corr(:,:,:,echo),vox);
    save_nii(nii,['src/mag_corr' num2str(echo) '.nii']);
    nii = make_nii(ph_corr(:,:,:,echo),vox);
    save_nii(nii,['src/ph_corr' num2str(echo) '.nii']);
end

% % do another BET on the mag_corr1.nii
% unix('bet2 src/mag_corr1.nii BET -f 0.4 -w 2 -m');
% unix('gunzip -f BET.nii.gz');
% unix('gunzip -f BET_mask.nii.gz');
% nii = load_nii('BET_mask.nii');
% mask = double(nii.img);


save('raw.mat','ph_corr','mag_corr','mask','-append');

% if need to clear variables for memory
clear mag_all ph_all


% (4) FUDGE laplacian unwrapping
disp('--> unwrap aliasing phase using fudge...');
for i = 1:imsize(4)
    unph(:,:,:,i) = fudge(ph_corr(:,:,:,i));
end
nii = make_nii(unph, vox);
save_nii(nii,'unph_fudge_before_correction.nii');


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

nii = make_nii(unph, vox);
save_nii(nii,'unph_fudge.nii');

unph_fudge = unph;
save('raw.mat','unph_fudge','-append');

mkdir('fudge');
cd('fudge');

clear unph_fudge ph_corr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set parameters
fit_thr = 20;
tik_reg = 1e-6;
cgs_num = 200;
tv_reg = 2e-4;
inv_num = 1000;


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

clear mag_corr

% extra filtering according to fitting residuals
% generate reliability map
fit_residual_0_blur = smooth3(fit_residual_0,'box',round(1./vox)*2+1); 
nii = make_nii(fit_residual_0_blur,vox);
save_nii(nii,'fit_residual_0_blur.nii');
R_0 = ones(size(fit_residual_0_blur));
R_0(fit_residual_0_blur >= fit_thr) = 0;

smv_rad = 3;
% RE-SHARP (tik_reg: Tikhonov regularization parameter)
disp('--> RESHARP to remove background field ...');
[lfs_resharp_0, mask_resharp_0] = resharp(tfs_0,mask.*R_0,vox,smv_rad,tik_reg,cgs_num);
% % 3D 2nd order polyfit to remove any residual background
% lfs_resharp= (lfs_resharp - poly3d(lfs_resharp,mask_resharp)).*mask_resharp;
% save nifti
mkdir('RESHARP');
nii = make_nii(lfs_resharp_0,vox);
save_nii(nii,['RESHARP/lfs_resharp_0_smvrad' num2str(smv_rad) '.nii']);

% iLSQR
chi_iLSQR_0 = QSM_iLSQR(lfs_resharp_0*(2.675e8*dicom_info.MagneticFieldStrength)/1e6,mask_resharp_0,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
nii = make_nii(chi_iLSQR_0,vox);
save_nii(nii,['RESHARP/chi_iLSQR_0_niter50_smvrad' num2str(smv_rad) '.nii']);



% perform inversion on each individual echo
tfs_resharp = zeros(imsize(1:4));
lfs_resharp = zeros(imsize(1:4));
chi_iLSQR = zeros(imsize(1:4));
mkdir('RESHARP');
for i = 1:imsize(4)
    % RE-SHARP (tik_reg: Tikhonov regularization parameter)
    disp('--> RESHARP to remove background field ...');
    [tfs_resharp(:,:,:,i), mask_resharp] = resharp(unph(:,:,:,i),mask,vox,smv_rad,tik_reg,cgs_num);
    lfs_resharp(:,:,:,i) = tfs_resharp(:,:,:,i)/TE(i)/(2.675e8*dicom_info.MagneticFieldStrength)*1e6; % unit ppm
    nii = make_nii(lfs_resharp(:,:,:,i),vox);
    save_nii(nii,['RESHARP/lfs_resharp_50' 'e' num2str(i) '.nii']);

    % iLSQR
    chi_iLSQR(:,:,:,i) = QSM_iLSQR(tfs_resharp(:,:,:,i),mask_resharp,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000*TE(i),'B0',dicom_info.MagneticFieldStrength);
    nii = make_nii(chi_iLSQR(:,:,:,i),vox);
    save_nii(nii,['RESHARP/chi_iLSQR_50' 'e' num2str(i) '.nii']);
end



save('cosmos.mat','lfs_resharp_0','mask_resharp_0','lfs_resharp','mask','vox','z_prjs');




%%%%%%%%%%%%%%%%%
% COSMOS
cd ../../../1/QSM_MEGE_7T/fudge/
load cosmos.mat
lfs_resharp1=lfs_resharp_0;
lfs_resharp_all1=lfs_resharp;
mask1=mask;
mask_resharp1=mask_resharp_0;
z_prjs1=z_prjs;
cd ../../../2/QSM_MEGE_7T/fudge/
load cosmos.mat
lfs_resharp2=lfs_resharp_0;
lfs_resharp_all2=lfs_resharp;
mask2=mask;
mask_resharp2=mask_resharp_0;
z_prjs2=z_prjs;
cd ../../../3/QSM_MEGE_7T/fudge/
load cosmos.mat
lfs_resharp3=lfs_resharp_0;
lfs_resharp_all3=lfs_resharp;
mask3=mask;
mask_resharp3=mask_resharp_0;
z_prjs3=z_prjs;
cd ../../../4/QSM_MEGE_7T/fudge/
load cosmos.mat
lfs_resharp4=lfs_resharp_0;
lfs_resharp_all4=lfs_resharp;
mask4=mask;
mask_resharp4=mask_resharp_0;
z_prjs4=z_prjs;
cd ../../../5/QSM_MEGE_7T/fudge/
load cosmos.mat
lfs_resharp5=lfs_resharp_0;
lfs_resharp_all5=lfs_resharp;
mask5=mask;
mask_resharp5=mask_resharp_0;
z_prjs5=z_prjs;
cd ../../../6/QSM_MEGE_7T/fudge/
load cosmos.mat
lfs_resharp6=lfs_resharp_0;
lfs_resharp_all6=lfs_resharp;
mask6=mask;
mask_resharp6=mask_resharp_0;
z_prjs6=z_prjs;



% define common space coordinates as the nature position orientation
% register all other orientations with the common space coordinate
/usr/local/fsl/5.0.9/bin/flirt -in /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/extension/QSM_MEGE_7T/RESHARP/chi_iLSQR_0_niter50_smvrad3.nii -ref /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/neutral/QSM_MEGE_7T/RESHARP/chi_iLSQR_0_niter50_smvrad3.nii -out /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/extension/flirt_qsm.nii.gz -omat /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/extension/flirt_qsm.mat -bins 256 -cost normcorr -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear

/usr/local/fsl/5.0.9/bin/flirt -in /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/flexion/QSM_MEGE_7T/RESHARP/chi_iLSQR_0_niter50_smvrad3.nii -ref /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/neutral/QSM_MEGE_7T/RESHARP/chi_iLSQR_0_niter50_smvrad3.nii -out /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/flexion/flirt_qsm.nii.gz -omat /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/flexion/flirt_qsm.mat -bins 256 -cost normcorr -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear

/usr/local/fsl/5.0.9/bin/flirt -in /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/right/QSM_MEGE_7T/RESHARP/chi_iLSQR_0_niter50_smvrad3.nii -ref /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/neutral/QSM_MEGE_7T/RESHARP/chi_iLSQR_0_niter50_smvrad3.nii -out /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/right/flirt_qsm.nii.gz -omat /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/right/flirt_qsm.mat -bins 256 -cost normcorr -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear

/usr/local/fsl/5.0.9/bin/flirt -in /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/left/QSM_MEGE_7T/RESHARP/chi_iLSQR_0_niter50_smvrad3.nii -ref /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/neutral/QSM_MEGE_7T/RESHARP/chi_iLSQR_0_niter50_smvrad3.nii -out /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/left/flirt_qsm.nii.gz -omat /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/left/flirt_qsm.mat -bins 256 -cost normcorr -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear



% apply the transformation to local field map
/usr/local/fsl/5.0.9/bin/flirt -in /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/extension/QSM_MEGE_7T/RESHARP/lfs_resharp_0_smvrad3_lsqr.nii -applyxfm -init /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/extension/flirt_qsm.mat -out /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/extension/flirt_lfs.nii -paddingsize 0.0 -interp trilinear -ref /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/neutral/QSM_MEGE_7T/RESHARP/lfs_resharp_0_smvrad3_lsqr.nii

/usr/local/fsl/5.0.9/bin/flirt -in /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/flexion/QSM_MEGE_7T/RESHARP/lfs_resharp_0_smvrad3_lsqr.nii -applyxfm -init /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/flexion/flirt_qsm.mat -out /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/flexion/flirt_lfs.nii -paddingsize 0.0 -interp trilinear -ref /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/neutral/QSM_MEGE_7T/RESHARP/lfs_resharp_0_smvrad3_lsqr.nii

/usr/local/fsl/5.0.9/bin/flirt -in /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/right/QSM_MEGE_7T/RESHARP/lfs_resharp_0_smvrad3_lsqr.nii -applyxfm -init /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/right/flirt_qsm.mat -out /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/right/flirt_lfs.nii -paddingsize 0.0 -interp trilinear -ref /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/neutral/QSM_MEGE_7T/RESHARP/lfs_resharp_0_smvrad3_lsqr.nii

/usr/local/fsl/5.0.9/bin/flirt -in /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/left/QSM_MEGE_7T/RESHARP/lfs_resharp_0_smvrad3_lsqr.nii -applyxfm -init /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/left/flirt_qsm.mat -out /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/left/flirt_lfs.nii -paddingsize 0.0 -interp trilinear -ref /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/neutral/QSM_MEGE_7T/RESHARP/lfs_resharp_0_smvrad3_lsqr.nii



% calculate the angles of B0 with registered local field maps
% R is transformation matrix from individual image space to common space (FLIRT matrix)
% common space coordinates = R* object image space coordinates
load /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/extension/flirt_qsm.mat -ASCII
R_t(:,:,1) = flirt_qsm(1:3,1:3);

load /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/flexion/flirt_qsm.mat -ASCII
R_t(:,:,2) = flirt_qsm(1:3,1:3);

load /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/right/flirt_qsm.mat -ASCII
R_t(:,:,3) = flirt_qsm(1:3,1:3);

load /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/left/flirt_qsm.mat -ASCII
R_t(:,:,4) = flirt_qsm(1:3,1:3);

R_t(:,:,5) = eye(3);


% (each orientation has own R and z_prjs)
% R is the rotation matrix from image space to common space
load('/gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/extension/QSM_MEGE_7T/raw.mat','z_prjs');
z_prjs_o(:,1) = z_prjs';

load('/gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/flexion/QSM_MEGE_7T/raw.mat','z_prjs');
z_prjs_o(:,2) = z_prjs';

load('/gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/right/QSM_MEGE_7T/raw.mat','z_prjs');
z_prjs_o(:,3) = z_prjs';

load('/gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/left/QSM_MEGE_7T/raw.mat','z_prjs');
z_prjs_o(:,4) = z_prjs';

load('/gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/neutral/QSM_MEGE_7T/raw.mat','z_prjs');
z_prjs_o(:,5) = z_prjs';

for i = 1:5
    z_prjs_c(:,i) = R_t(:,:,i)'*z_prjs_o(:,i);
end


%% COSMOS reconstruction with closed-form solution
% load in registered local field shift maps
unix('gunzip -f /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/extension/*.gz');
nii = load_nii('/gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/extension/flirt_lfs.nii');
lfs(:,:,:,1) = double(nii.img);

unix('gunzip -f /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/flexion/*.gz');
nii = load_nii('/gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/flexion/flirt_lfs.nii');
lfs(:,:,:,2) = double(nii.img);

unix('gunzip -f /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/right/*.gz');
nii = load_nii('/gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/right/flirt_lfs.nii');
lfs(:,:,:,3) = double(nii.img);

unix('gunzip -f /gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/left/*.gz');
nii = load_nii('/gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/left/flirt_lfs.nii');
lfs(:,:,:,4) = double(nii.img);

nii = load_nii('/gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/neutral/QSM_MEGE_7T/RESHARP/lfs_resharp_0_smvrad3_lsqr.nii');
lfs(:,:,:,5) = double(nii.img);

mask = and(and(and(and(lfs(:,:,:,1),lfs(:,:,:,2)),lfs(:,:,:,3)),lfs(:,:,:,4)),lfs(:,:,:,5));
mask = double(mask);

% construct k-space kernel for each orientation
% create K-space filter kernel D
%%%%% make this a seperate function in the future
load('/gpfs/M2Scratch/NCIgb5/hongfu/recon/ED/neutral/QSM_MEGE_7T/raw.mat','vox');

Nx = size(lfs,1);
Ny = size(lfs,2);
Nz = size(lfs,3);

for i = 1:5
FOV = vox.*[Nx,Ny,Nz];
FOVx = FOV(1);
FOVy = FOV(2);
FOVz = FOV(3);

x = -Nx/2:Nx/2-1;
y = -Ny/2:Ny/2-1;
z = -Nz/2:Nz/2-1;
[kx,ky,kz] = ndgrid(x/FOVx,y/FOVy,z/FOVz);
% D = 1/3 - kz.^2./(kx.^2 + ky.^2 + kz.^2);
D = 1/3 - (kx.*z_prjs_c(1,i)+ky.*z_prjs_c(2,i)+kz.*z_prjs_c(3,i)).^2./(kx.^2 + ky.^2 + kz.^2);
D(floor(Nx/2+1),floor(Ny/2+1),floor(Nz/2+1)) = 0;
D = fftshift(D);

kernel(:,:,:,i) = D;
end


for i = 1:5
    lfs_k(:,:,:,i) = fftn(lfs(:,:,:,i).*mask);
end

kernel_sum = sum(abs(kernel).^2, 4);

chi_cosmos = real( ifftn( sum(kernel .* lfs_k, 4) ./ (eps + kernel_sum) ) ) .* mask;

nii = make_nii(chi_cosmos,vox);
save_nii(nii,'cosmos_5.nii');

