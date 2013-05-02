function sus = qsm(path_in, path_out, params)
%QSM Quantitative susceptibility mapping.
%   SUS = QSM(PATH_IN, PATH_OUT, PARAMS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH_IN    - directory of .fid from gemsme3d sequence  : pwd/gemsme3d_R2s_01.fid
%   PATH_OUT   - directory to save nifti and/or matrixes   : pwd
%   PARAMS     - parameter structure including fields below (!in small case!)
%    .ker_rad  - radius (mm) of RESHARP convolution kernel : 5
%    .tik_reg  - Tikhonov regularization parameter         : 0.005
%    .tv_reg   - Total variation regularization parameter  : 0.001
%    .save_mat - whether to save matrixes (1) or not (0)   : 1
%   SUS        - susceptibility maps as output

%% default settings and prompts
if ~ exist('path_in','var') || isempty(path_in)
    path_in = [pwd '/gemsme3d_R2s_01.fid'];
end

if ~ exist('path_out','var') || isempty(path_out)
    path_out = pwd;
end

if ~ exist('params','var') || isempty(params)
    params = [];
end

if ~ isfield(params,'ker_rad')
    params.ker_rad = 5;
end

if ~ isfield(params,'tik_reg')
    params.tik_reg = 5e-3;
end

if ~ isfield(params,'tv_reg')
    params.tv_reg = 1e-3;
end

if ~ isfield(params,'save_mat')
    params.save_mat = 1;
end

ker_rad  = params.ker_rad;
tik_reg  = params.tik_reg;
tv_reg   = params.tv_reg;
save_mat = params.save_mat;

fprintf(['\n' ...
'Please confirm the settings of the following required parameters \n\n' ...
'--> path_in  -- directory of gemsme3d_R2s_01.fid rawdata  :  %s \n' ...
'--> path_out -- directory to save nifti and/or matrixes   :  %s \n' ...
'--> ker_rad  -- radius (mm) of RESHARP convolution kernel :  %g \n' ...
'--> tik_reg  -- Tikhonov regularization parameter         :  %g \n' ...
'--> tv_reg   -- Total variation regularization parameter  :  %g \n' ...
'--> save_mat -- whether to save matrixes (1) or not (0)   :  %g \n\n' ...
'Start in 10 sec, Ctrl-C to terminate!\n\n'], ...
path_in, path_out, ker_rad, tik_reg, tv_reg, save_mat);

pause(10);


%% define directories
if save_mat
    path_mat = [path_out '/matrix'];
    mkdir(path_mat);
end

path_nft = [path_out '/nifti'];
mkdir(path_nft);

TEMP = [path_out '/temp'];
mkdir(TEMP);

cd(TEMP);


%% reconstruct complex image from fid file
disp('--> (1/9) reconstruct fid to complex img ...');
[img,par] = reconfid(path_in);

% keep only the first 5 echoes
img = img(:,:,:,1:5,:);
par.ne = 5;
[np,nv,nv2,ne,~] = size(img);
res = par.res; % resolution in mm/pixel


% save matrix
if save_mat
    mkdir([path_mat '/rawdata']);
    save([path_mat '/rawdata/img.mat'],'img','-v7.3');
    save([path_mat '/rawdata/par.mat'],'par','-v7.3');
end


%% combine magnitude channels using marc's 'arrayrec.m'
disp('--> (2/9) combine rcvrs for magnitude ...');
mag_cmb = zeros(np,nv,nv2,ne);
for echo = 1:ne
    mag_cmb(:,:,:,echo) = arrayrec(squeeze(img(:,:,:,echo,:)),1/2);
end

% save matrix
if save_mat
    mkdir([path_mat '/combine']);
    save([path_mat '/combine/mag_cmb.mat'],'mag_cmb','-v7.3');
end

% save nifti
mkdir([path_nft '/combine']);
for echo =  1:ne
    nii = make_nii(mag_cmb(:,:,:,echo),res);
    save_nii(nii,[path_nft '/combine/mag_te' num2str(echo) '.nii']);
end


%% generate mask from SOS combined magnitude of first echo
disp('--> (3/9) extract brain volume and generate mask ...');
setenv('path_nft', path_nft);
unix('bet $path_nft/combine/mag_te1.nii BET -f 0.5 -m -R');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);

% save matrix
if save_mat
    mkdir([path_mat '/mask']);
    save([path_mat '/mask/mask.mat'],'mask','-v7.3');
end

% save nifti
mkdir([path_nft '/mask']);
copyfile('BET_mask.nii',[path_nft '/mask/mask.nii']);


%% combine phase channels
disp('--> (4/9) combine rcvrs for phase ...');
ph_cmb = sense(img,par);

% save matrix
if save_mat
    save([path_mat '/combine/ph_cmb.mat'],'ph_cmb','-v7.3');
end

% save nifti
for echo =  1:ne
    nii = make_nii(ph_cmb(:,:,:,echo),res);
    save_nii(nii,[path_nft '/combine/ph_te' num2str(echo) '.nii']);
end

clear img;


%% unwrap phase from each echo
disp('--> (5/9) unwrap aliasing phase for each TE ...');

bash_command = sprintf(['for ph in $path_nft/combine/ph*\n' ...
'do\n' ...
'	base=`basename $ph`;\n' ...
'	dir=`dirname $ph`;\n' ...
'	mag=$dir/"mag"${base:2};\n' ...
'	unph="unph"${base:2};\n' ...
'	prelude -a $mag -p $ph -u $unph -m BET_mask.nii -n 8&\n' ...
'done\n' ...
'wait\n' ...
'gunzip -f unph*.gz\n']);

unix(bash_command);

unph_cmb = zeros(np,nv,nv2,ne);
for echo = 1:ne
    nii = load_nii(['unph_te' num2str(echo) '.nii']);
    unph_cmb(:,:,:,echo) = double(nii.img);
end


% check and correct for 2pi jump between echoes
disp('--> (6/9) correct for potential 2pi jumps between TEs ...')

nii = load_nii('unph_diff.nii');
unph_diff = double(nii.img);

for echo = 2:ne
    meandiff = unph_cmb(:,:,:,echo)-unph_cmb(:,:,:,1)-(echo-1)*unph_diff;
    meandiff = meandiff(mask==1);
    meandiff = mean(meandiff(:));
    njump = round(meandiff/(2*pi));
    disp(['    ' num2str(njump) ' 2pi jumps for TE' num2str(echo)]);
    unph_cmb(:,:,:,echo) = unph_cmb(:,:,:,echo) - njump*2*pi;
end

% save matrix
if save_mat
    mkdir([path_mat '/unwrap']);
    save([path_mat '/unwrap/unph_cmb.mat'],'unph_cmb','-v7.3'); 
end

% save nifti
mkdir([path_nft '/unwrap']);
for echo = 1:ne
    nii = make_nii(unph_cmb(:,:,:,echo),res);
    save_nii(nii,[path_nft '/unwrap/unph_te' num2str(echo) '.nii']);
end

clear ph_cmb


%% fit phase images with echo times
disp('--> (7/9) magnitude weighted LS fit of phase to TE ...');
[tfs, fit_residual] = echofit(unph_cmb,mag_cmb,par); 

% generate reliability map
R = ones(size(fit_residual));
R(fit_residual >= 10) = 0;

% save matrix
if save_mat
    mkdir([path_mat '/fit']);
    save([path_mat '/fit/tfs.mat'],'tfs','-v7.3');
    save([path_mat '/fit/fit_residual.mat'],'fit_residual','-v7.3');
    save([path_mat '/fit/R.mat'],'R','-v7.3');
end

% save nifti
mkdir([path_nft '/fit']);
nii = make_nii(tfs,res);
save_nii(nii,[path_nft '/fit/tfs.nii']);
nii = make_nii(fit_residual,res);
save_nii(nii,[path_nft '/fit/fit_residual.nii']);
nii = make_nii(R,res);
save_nii(nii,[path_nft '/fit/R.nii']);

clear unph_cmb


%% RESHARP (tik_reg: Tikhonov regularization parameter)
disp('--> (8/9) RESHARP to remove background field ...');
[lfs, mask_ero] = resharp(tfs,mask.*R,par,ker_rad,tik_reg);

% save matrix
if save_mat
    mkdir([path_mat '/rmbkg']);
    save([path_mat '/rmbkg/lfs.mat'],'lfs','-v7.3'); 
    save([path_mat '/mask/mask_ero.mat'],'mask_ero','-v7.3'); 
end

% save nifti
mkdir([path_nft '/rmbkg/']);
nii = make_nii(lfs,res);
save_nii(nii,[path_nft '/rmbkg/lfs_xy_' num2str(tik_reg) '.nii']);
nii = make_nii(permute(lfs,[1 3 2]),res);
save_nii(nii,[path_nft '/rmbkg/lfs_xz_' num2str(tik_reg) '.nii']);
nii = make_nii(permute(lfs,[2 3 1]),res);
save_nii(nii,[path_nft '/rmbkg/lfs_yz_' num2str(tik_reg) '.nii']);

nii = make_nii(mask_ero,res);
save_nii(nii,[path_nft '/mask/mask_ero.nii']);


%% inversion of RESHARP
disp('--> (9/9) Total variation susceptibility inversion ...');
sus = tvdi(lfs, mask_ero, par, tv_reg, mag_cmb(:,:,:,5)); 

% save matrix
if save_mat
    mkdir([path_mat '/inversion']);
    save([path_mat '/inversion/sus.mat'],'sus','-v7.3'); 
end

% save nifti
mkdir([path_nft '/inversion']);
nii = make_nii(sus,res);
save_nii(nii,[path_nft '/inversion/sus_xy_' num2str(tv_reg) '.nii']);
nii = make_nii(permute(sus,[1 3 2]),res);
save_nii(nii,[path_nft '/inversion/sus_xz_' num2str(tv_reg) '.nii']);
nii = make_nii(permute(sus,[2 3 1]),res);
save_nii(nii,[path_nft '/inversion/sus_yz_' num2str(tv_reg) '.nii']);


