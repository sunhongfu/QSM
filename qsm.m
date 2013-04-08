%QSM Quantitative susceptibility mapping.
%   qsm
%
%   PATH_IN  : directory of .fid from 'gemsme3d' sequence
%   PATH_OUT : directory for outputs
%   KER_RAD  : the radius of the convolution kernel (mm)
%   TIK_REG  : Tikhonov regularization parameter for RESHARP
%   TV_REG   : Total variation regularization parameter for inversion
%   SAVE_MAT : to save all the matrixes


%% default settings
if ~ exist('PATH_IN','var')
    PATH_IN = pwd;
end

if ~ exist('PATH_OUT','var')
    PATH_OUT = pwd;
end

if ~ exist('KER_RAD','var')
    KER_RAD = 5;
end

if ~ exist('TIK_REG','var')
    TIK_REG = 5e-3;
end

if ~ exist('TV_REG','var')
    TV_REG = 5e-4;
end

if ~ exist('SAVE_MAT','var')
    SAVE_MAT = 1;
end

% comment this following part if you do not need confirmation prompts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
confirm = input([ ...
'***Please confirm the settings of the following required parameters***\n\n' ...
'(1) PATH_IN  - directory of .fid from gemsme3d sequence  : ' PATH_IN '\n' ...
'(2) PATH_OUT - directory to save nifti and/or matrixes   : ' PATH_OUT '\n' ...
'(3) KER_RAD  - radius (mm) of RESHARP convolution kernel : ' num2str(KER_RAD) '\n' ...
'(4) TIK_REG  - Tikhonov regularization parameter         : ' num2str(TIK_REG) '\n' ...
'(5) TV_REG   - Total variation regularization parameter  : ' num2str(TV_REG) '\n' ...
'(6) SAVE_MAT - whether to save matrixes (1) or not (0)   : ' num2str(SAVE_MAT) '\n\n' ...
'To confirm the above parameters, enter "y", otherwise "n": '], 's');

while true
    if lower(confirm) == 'y'
        disp('Begin QSM reconstruction!');
        break;
    elseif lower(confirm) == 'n'
        error('Please re-define the above parameters and re-run the srcipt');
    else
        confirm = input('please enter "y" for yes, or "n" for no: ', 's');
        continue;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% define directories
MATRIX = [PATH_OUT '/matrix'];
NIFTI  = [PATH_OUT '/nifti'];
TMP    = [PATH_OUT '/tmp'];

mkdir(MATRIX);
mkdir(NIFTI);
mkdir(TMP);

cd(TMP);


%% reconstruct complex image from fid file
disp('--> (1/9) recon fid to complex img ...');
[img,par] = reconfid(PATH_IN);
[np,nv,nv2,ne,~] = size(img);
res = par.res; % resolution in mm/pixel

% save matrix
if SAVE_MAT
    mkdir([MATRIX '/rawdata']);
    save([MATRIX '/rawdata/img.mat'],'img','-v7.3');
    save([MATRIX '/rawdata/par.mat'],'par','-v7.3');
end


%% combine magnitude channels using marc's 'arrayrec.m'
disp('--> (2/9) combine rcvrs for magnitude ...');
mag_cmb = zeros(np,nv,nv2,ne);
for echo = 1:ne
    mag_cmb(:,:,:,echo) = arrayrec(squeeze(img(:,:,:,echo,:)),1/2);
end

% save matrix
if SAVE_MAT
    mkdir([MATRIX '/cmb']);
    save([MATRIX '/cmb/mag_cmb.mat'],'mag_cmb','-v7.3');
end

% save nifti
mkdir([NIFTI '/cmb']);
for echo =  1:ne
    nii = make_nii(mag_cmb(:,:,:,echo),res);
    save_nii(nii,[NIFTI '/cmb/mag_te' num2str(echo) '.nii']);
end


%% generate mask from SOS combined magnitude of first echo
disp('--> (3/9) extract brain volume and generate mask ...');
setenv('NIFTI', NIFTI);
unix('bet $NIFTI/cmb/mag_te1.nii BET -f 0.5 -m -R');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);

% save matrix
if SAVE_MAT
    mkdir([MATRIX '/mask']);
    save([MATRIX '/mask/mask.mat'],'mask','-v7.3');
end

% save nifti
mkdir([NIFTI '/mask']);
copyfile('BET_mask.nii',[NIFTI '/mask/mask.nii']);


%% combine phase channels using MCPC-3D
disp('--> (4/9) combine rcvrs for phase ...');
ph_cmb = mcpc3d(img,par);

% save matrix
if SAVE_MAT
    save([MATRIX '/cmb/ph_cmb.mat'],'ph_cmb','-v7.3');
end

% save nifti
for echo =  1:ne
    nii = make_nii(ph_cmb(:,:,:,echo),res);
    save_nii(nii,[NIFTI '/cmb/ph_te' num2str(echo) '.nii']);
end

clear img;


%% unwrap phase from each echo
disp('--> (5/9) initial unwrap aliasing phase for each TE ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unix('~/projects/QSM/scripts/unwrap $NIFTI/cmb/ph*');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unph = zeros(np,nv,nv2,ne);
for echo = 1:ne
    nii = load_nii(['unph_te' num2str(echo) '.nii']);
    unph(:,:,:,echo) = double(nii.img);
end

% check and correct for 2pi jump between echoes
disp('--> (6/9) correct for potential 2pi jumps between TEs ...')
%****************************************%
% might be better to unwrap the phase_diff
unph_diff = unph(:,:,:,2) - unph(:,:,:,1);
%****************************************%
for echo = 3:ne
    meandiff = unph(:,:,:,echo)-unph(:,:,:,1)-(echo-1)*unph_diff;
    meandiff = meandiff(mask==1);
    meandiff = median(meandiff(:));
    njump = round(meandiff/(2*pi));
    disp([num2str(njump) ' 2pi jumps for TE' num2str(echo)]);
    unph(:,:,:,echo) = unph(:,:,:,echo) - njump*2*pi;
end

% save matrix
if SAVE_MAT
    mkdir([MATRIX '/unwrap']);
    save([MATRIX '/unwrap/unph.mat'],'unph','-v7.3'); 
end

% save nifti
mkdir([NIFTI '/unwrap']);
for echo = 1:ne
    nii = make_nii(unph(:,:,:,echo),res);
    save_nii(nii,[NIFTI '/unwrap/unph_te' num2str(echo) '.nii']);
end

clear ph_cmb unph_diff


%% fit phase images with echo times
disp('--> (7/9) magnitude weighted LS fit of phase to TE ...');
[field_total, fit_residual] = echofit(unph(:,:,:,1:5),mag_cmb(:,:,:,1:5),par); 

% generate reliability map
R = ones(size(fit_residual));
R(fit_residual >= 1) = 0;

% save matrix
if SAVE_MAT
    mkdir([MATRIX '/fit']);
    save([MATRIX '/fit/field_total.mat'],'field_total','-v7.3');
    save([MATRIX '/fit/fit_residual.mat'],'fit_residual','-v7.3');
    save([MATRIX '/fit/R.mat'],'R','-v7.3');
end

% save nifti
mkdir([NIFTI '/fit']);
nii = make_nii(field_total,res);
save_nii(nii,[NIFTI '/fit/field_total.nii']);
nii = make_nii(fit_residual,res);
save_nii(nii,[NIFTI '/fit/fit_residual.nii']);
nii = make_nii(R,res);
save_nii(nii,[NIFTI '/fit/R.nii']);

clear unph mag_cmb fit_residual


%% RESHARP (TIK_REG: Tikhonov regularization parameter)
disp('--> (8/9) RESHARP to remove background field ...');
[field_local, mask_ero] = resharp(field_total.*R,mask.*R,par,KER_RAD,TIK_REG);

% save matrix
if SAVE_MAT
    mkdir([MATRIX '/rmbkg/resharp']);
    save([MATRIX '/rmbkg/resharp/field_local.mat'],'field_local','-v7.3'); 
    save([MATRIX '/rmbkg/resharp/mask_ero.mat'],'mask_ero','-v7.3'); 
end

% save nifti
mkdir([NIFTI '/rmbkg/resharp/']);
nii = make_nii(field_local,res);
save_nii(nii,[NIFTI '/rmbkg/resharp/resharp_' num2str(TIK_REG) '.nii']);
nii = make_nii(mask_ero,res);
save_nii(nii,[NIFTI '/rmbkg/resharp/mask_ero.nii']);

clear field_total mask R


%% inversion of SHARP
disp('--> (9/9) Total variation susceptibility inversion ...');
sus = tvdi(field_local, mask_ero, par, TV_REG); 

% save matrix
if SAVE_MAT
    mkdir([MATRIX '/inversion']);
    save([MATRIX '/inversion/sus.mat'],'sus','-v7.3'); 
end

% save nifti
mkdir([NIFTI '/inversion']);
nii = make_nii(sus,res);
save_nii(nii,[NIFTI '/inversion/sus_xy_' num2str(TV_REG) '.nii']);
nii = make_nii(permute(sus,[1 3 2]),res);
save_nii(nii,[NIFTI '/inversion/sus_xz_' num2str(TV_REG) '.nii']);
nii = make_nii(permute(sus,[2 3 1]),res);
save_nii(nii,[NIFTI '/inversion/sus_yz_' num2str(TV_REG) '.nii']);


