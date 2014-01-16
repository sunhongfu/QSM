function qsm_swi47(path_in, path_out, options)
%QSM_SWI47 Quantitative susceptibility mapping from SWI sequence at 4.7T.
%   QSM_SWI47(PATH_IN, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH_IN    - directory of .fid from ge3d sequence      : ge3d__01.fid
%   PATH_OUT   - directory to save nifti and/or matrixes   : QSM_SWI_vxxx
%   OPTIONS    - parameter structure including fields below
%    .ref_coi  - reference coil to use for phase combine   : 3
%    .eig_rad  - radius (mm) of eig decomp kernel          : 3
%    .smv_rad  - radius (mm) of SMV convolution kernel     : 6
%    .tik_reg  - Tikhonov regularization for RESHARP       : 0.001
%    .tv_reg   - Total variation regularization parameter  : 0.0005
%    .bet_thr  - threshold for BET brain mask              : 0.4
%    .tvdi_n   - iteration number of TVDI (nlcg)           : 200


%% default settings
if ~ exist('path_in','var') || isempty(path_in)
    path_in = pwd;
end

if exist([path_in '/fid'],'file')
    path_fid = path_in;
elseif exist([path_in '/ge3d__01.fid/fid'],'file')
    path_fid = [path_in, '/ge3d__01.fid'];
else
    error('cannot find .fid file');
end

if ~ exist('path_out','var') || isempty(path_out)
    path_out = path_fid;
end

if ~ exist('options','var') || isempty(options)
    options = [];
end

if ~ isfield(options,'ref_coi')
    options.ref_coi = 3;
end

if ~ isfield(options,'eig_rad')
    options.eig_rad = 3;
end

if ~ isfield(options,'bet_thr')
    options.bet_thr = 0.4;
end

if ~ isfield(options,'smv_rad')
    options.smv_rad = 6;
end

if ~ isfield(options,'tik_reg')
    options.tik_reg = 1e-3;
end

if ~ isfield(options,'tv_reg')
    options.tv_reg = 5e-4;
end

if ~ isfield(options,'tvdi_n')
    options.tvdi_n = 200;
end

ref_coi = options.ref_coi;
eig_rad = options.eig_rad;
bet_thr = options.bet_thr;
smv_rad = options.smv_rad;
tik_reg = options.tik_reg;
tv_reg  = options.tv_reg;
tvdi_n  = options.tvdi_n;


%% define directories
path_qsm = [path_out '/QSM_SWI_v200'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);


%% generate raw img
disp('--> reconstruct fid to complex img ...');
[img,Pars] = swi47_recon(path_fid);


%% interpolate to iso-resoluation in plane
k = ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(img,1),2),[],1),[],2),1),2);
pad = round((Pars.np/2 * Pars.lpe / Pars.lro - Pars.nv)/2);
k = padarray(k,[0 pad]);
img = fftshift(fftshift(fft(fft(fftshift(fftshift(k,1),2),[],1),[],2),1),2);

% swap the first two dimensions
img = permute(img, [2 1 3 4]);
[nv,np,ns,~] = size(img); % phase, readout, slice, receivers
vox = [Pars.lpe/nv, Pars.lro/np, Pars.lpe2/ns]*10;


%% combine RF coils
disp('--> combine RF rcvrs ...');
img_cmb = sense_se(img,vox,ref_coi,eig_rad);

% save nifti
mkdir('combine');
nii = make_nii(abs(img_cmb),vox);
save_nii(nii,'combine/mag_cmb.nii');

clear img;


%% center k-space correction (readout direction)
k = ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(img_cmb,1),2),[],1),[],2),1),2);
[~,Ind] = max(abs(k(:)));
Ix = ceil(mod(Ind,np*nv)/nv);

% Apply phase ramp
pix = np/2-Ix; % voxel shift
ph_ramp = exp(-sqrt(-1)*2*pi*pix*(-1/2:1/np:1/2-1/np));
img_cmb = img_cmb.* repmat(ph_ramp,[nv 1 ns]);

% save nifti
nii = make_nii(angle(img_cmb),vox);
save_nii(nii,'combine/ph_cmb.nii');


%% generate brain mask
disp('--> extract brain volume and generate mask ...');
setenv('path_qsm', path_qsm);
setenv('bet_thr',num2str(bet_thr));
unix('bet $path_qsm/combine/mag_cmb.nii BET -f ${bet_thr} -m -R');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);


%% unwrap combined phase with PRELUDE
disp('--> unwrap aliasing phase ...');
unix('prelude -a $path_qsm/combine/mag_cmb.nii -p $path_qsm/combine/ph_cmb.nii -u unph.nii -m BET_mask.nii -n 8');
unix('gunzip -f unph.nii.gz');
nii = load_nii('unph.nii');
unph = double(nii.img);


%% background field removal
disp('--> RESHARP to remove background field ...');
mkdir('RESHARP');
[lph_resharp,mask_resharp] = resharp(unph,mask,vox,smv_rad,tik_reg);

% normalize to echo time and field strength
% ph = gamma*dB*TE
% dB/B = ph/(gamma*TE*B0)
% units: TE s, gamma 2.675e8 rad/(sT), B0 4.7T
lfs_resharp = -lph_resharp/(2.675e8*Pars.te*4.7)*1e6; % unit ppm

nii = make_nii(lfs_resharp,vox);
save_nii(nii,'RESHARP/lfs_resharp.nii');


%% susceptibility inversion
disp('--> TV susceptibility inversion ...');
sus_resharp = tvdi(lfs_resharp,mask_resharp,vox,tv_reg,abs(img_cmb),tvdi_n);

nii = make_nii(sus_resharp,vox);
save_nii(nii,'RESHARP/sus_resharp.nii');


%% clean up
clear nii;
save('all.mat','-v7.3');
unix('rm *.nii*');
cd(init_dir);

