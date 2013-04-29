%QSM Quantitative susceptibility mapping.
%   QSM is the main script to reconstruct QSM from the R2* sequence.
%
%   The following parameter settings need to be re-defined if necessary:
%   (1) PATH_IN  - directory of .fid from gemsme3d sequence  : pwd
%   (2) PATH_OUT - directory to save nifti and/or matrixes   : pwd
%   (3) KER_RAD  - radius (mm) of RESHARP convolution kernel : 5
%   (4) TIK_REG  - Tikhonov regularization parameter         : 0.005
%   (5) TV_REG   - Total variation regularization parameter  : 0.0005
%   (6) SAVE_MAT - whether to save matrixes (1) or not (0)   : 1


%% default settings and prompts
if ~ exist('PATH_IN','var')
    PATH_IN = [pwd '/gemsme3d_R2s_01.fid'];
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
    TV_REG = 1e-3;
end

if ~ exist('SAVE_MAT','var')
    SAVE_MAT = 1;
end

% comment this following part if you do not need confirmation prompts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
confirm = input([ ...
'*** Please confirm the settings of the following required parameters ***\n\n' ...
'(1) PATH_IN  - directory of gemsme3d_R2s_01.fid rawdata  : ' PATH_IN '\n' ...
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
        error('Please re-define any above parameters and re-run the srcipt');
    else
        confirm = input('please enter "y" for yes, or "n" for no: ', 's');
        continue;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% define directories
if SAVE_MAT
    MATRIX = [PATH_OUT '/matrix'];
    mkdir(MATRIX);
end

NIFTI = [PATH_OUT '/nifti'];
mkdir(NIFTI);

TEMP = [PATH_OUT '/temp'];
mkdir(TEMP);

cd(TEMP);


%% reconstruct complex image from fid file
disp('--> (1/9) reconstruct fid to complex img ...');
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
    mkdir([MATRIX '/combine']);
    save([MATRIX '/combine/mag_cmb.mat'],'mag_cmb','-v7.3');
end

% save nifti
mkdir([NIFTI '/combine']);
for echo =  1:ne
    nii = make_nii(mag_cmb(:,:,:,echo),res);
    save_nii(nii,[NIFTI '/combine/mag_te' num2str(echo) '.nii']);
end


%% generate mask from SOS combined magnitude of first echo
disp('--> (3/9) extract brain volume and generate mask ...');
setenv('NIFTI', NIFTI);
unix('bet $NIFTI/combine/mag_te1.nii BET -f 0.5 -m -R');
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
    save([MATRIX '/combine/ph_cmb.mat'],'ph_cmb','-v7.3');
end

% save nifti
for echo =  1:ne
    nii = make_nii(ph_cmb(:,:,:,echo),res);
    save_nii(nii,[NIFTI '/combine/ph_te' num2str(echo) '.nii']);
end

clear img;


%% unwrap phase from each echo
disp('--> (5/9) unwrap aliasing phase for each TE ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% external bash script
unix('~/projects/QSM/scripts/unwrap $NIFTI/combine/ph*');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unph_cmb = zeros(np,nv,nv2,ne);
for echo = 1:ne
    nii = load_nii(['unph_te' num2str(echo) '.nii']);
    unph_cmb(:,:,:,echo) = double(nii.img);
end

% check and correct for 2pi jump between echoes
disp('--> (6/9) correct for potential 2pi jumps between TEs ...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) calucate the unph_diff directly (uncomment the following line)
unph_diff = unph_cmb(:,:,:,2) - unph_cmb(:,:,:,1);
% (2) might be better to unwrap the ph_diff
% ph_diff = angle(exp(1j*(unph_cmb(:,:,:,2) - unph_cmb(:,:,:,1))));
% nii = make_nii(ph_diff,res);
% save_nii(nii,'ph_diff.nii');
% unix('prelude -a BET.nii -p ph_diff.nii -u unph_diff.nii -m BET_mask.nii -n 8');
% unix('gunzip -f unph_diff.nii.gz');
% nii = load_nii('unph_diff.nii');
% unph_diff = double(nii.img);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for echo = 2:ne
    meandiff = unph_cmb(:,:,:,echo)-unph_cmb(:,:,:,1)-(echo-1)*unph_diff;
    meandiff = meandiff(mask==1);
    meandiff = mean(meandiff(:));
    njump = round(meandiff/(2*pi));
    disp(['    ' num2str(njump) ' 2pi jumps for TE' num2str(echo)]);
    unph_cmb(:,:,:,echo) = unph_cmb(:,:,:,echo) - njump*2*pi;
end

% save matrix
if SAVE_MAT
    mkdir([MATRIX '/unwrap']);
    save([MATRIX '/unwrap/unph_cmb.mat'],'unph_cmb','-v7.3'); 
end

% save nifti
mkdir([NIFTI '/unwrap']);
for echo = 1:ne
    nii = make_nii(unph_cmb(:,:,:,echo),res);
    save_nii(nii,[NIFTI '/unwrap/unph_te' num2str(echo) '.nii']);
end

clear ph_cmb


%% fit phase images with echo times
disp('--> (7/9) magnitude weighted LS fit of phase to TE ...');
[tfs, fit_residual] = echofit(unph_cmb(:,:,:,1:5),mag_cmb(:,:,:,1:5),par); 

% generate reliability map
R = ones(size(fit_residual));
R(fit_residual >= 10) = 0;

% save matrix
if SAVE_MAT
    mkdir([MATRIX '/fit']);
    save([MATRIX '/fit/tfs.mat'],'tfs','-v7.3');
    save([MATRIX '/fit/fit_residual.mat'],'fit_residual','-v7.3');
    save([MATRIX '/fit/R.mat'],'R','-v7.3');
end

% save nifti
mkdir([NIFTI '/fit']);
nii = make_nii(tfs,res);
save_nii(nii,[NIFTI '/fit/tfs.nii']);
nii = make_nii(fit_residual,res);
save_nii(nii,[NIFTI '/fit/fit_residual.nii']);
nii = make_nii(R,res);
save_nii(nii,[NIFTI '/fit/R.nii']);

clear unph_cmb mag_cmb


%% RESHARP (TIK_REG: Tikhonov regularization parameter)
disp('--> (8/9) RESHARP to remove background field ...');
[lfs, mask_ero] = resharp(tfs,mask.*R,par,KER_RAD,TIK_REG);

% save matrix
if SAVE_MAT
    mkdir([MATRIX '/rmbkg']);
    save([MATRIX '/rmbkg/lfs.mat'],'lfs','-v7.3'); 
    save([MATRIX '/mask/mask_ero.mat'],'mask_ero','-v7.3'); 
end

% save nifti
mkdir([NIFTI '/rmbkg/']);
nii = make_nii(lfs,res);
save_nii(nii,[NIFTI '/rmbkg/lfs_xy_' num2str(TIK_REG) '.nii']);
nii = make_nii(permute(lfs,[1 3 2]),res);
save_nii(nii,[NIFTI '/rmbkg/lfs_xz_' num2str(TIK_REG) '.nii']);
nii = make_nii(permute(lfs,[2 3 1]),res);
save_nii(nii,[NIFTI '/rmbkg/lfs_yz_' num2str(TIK_REG) '.nii']);

nii = make_nii(mask_ero,res);
save_nii(nii,[NIFTI '/mask/mask_ero.nii']);


%% inversion of RESHARP
disp('--> (9/9) Total variation susceptibility inversion ...');
sus = tvdi(lfs, mask_ero, par, TV_REG); 

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


