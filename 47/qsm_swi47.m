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
%    .smv_rad  - radius (mm) of SMV convolution kernel     : 4
%    .tik_reg  - Tikhonov regularization for resharp       : 0.0005
%    .tv_reg   - Total variation regularization parameter  : 0.0005
%    .bet_thr  - threshold for BET brain mask              : 0.3
%    .tvdi_n   - iteration number of TVDI (nlcg)           : 200
%    .sav_all  - save all the variables for debug          : 1


%%% default settings
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
    options.bet_thr = 0.3;
end

if ~ isfield(options,'smv_rad')
    options.smv_rad = 4;
end

if ~ isfield(options,'tik_reg')
    options.tik_reg = 5e-4;
end

if ~ isfield(options,'tv_reg')
    options.tv_reg = 5e-4;
end

if ~ isfield(options,'tvdi_n')
    options.tvdi_n = 200;
end

if ~ isfield(options,'sav_all')
    options.sav_all = 1;
end

ref_coi = options.ref_coi;
eig_rad = options.eig_rad;
bet_thr = options.bet_thr;
smv_rad = options.smv_rad;
tik_reg = options.tik_reg;
tv_reg  = options.tv_reg;
tvdi_n  = options.tvdi_n;
sav_all = options.sav_all;



%%% define directories
path_qsm = [path_out '/QSM_SWI_v300'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);


%%% generate raw img
disp('--> reconstruct fid to complex img ...');
[img,Pars] = swi47_recon(path_fid);


%%% interpolate to iso-resoluation in plane
% k = ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(img,1),2),[],1),[],2),1),2);
% pad = round((Pars.np/2 * Pars.lpe / Pars.lro - Pars.nv)/2);
% k = padarray(k,[0 pad]);
% img = fftshift(fftshift(fft(fft(fftshift(fftshift(k,1),2),[],1),[],2),1),2);

k = fft(fft(img,[],1),[],2);
pad = round(Pars.np/2*Pars.lpe / Pars.lro - Pars.nv);
imsize = size(k);
if mod(imsize(2),2) % if size of k is odd
    k_pad = ifftshift(padarray(padarray(fftshift(k,2),[0 round(pad/2)],'pre'), [0 pad-round(pad/2)], 'post'),2);
else % size of k is even
    k_s = fftshift(k,2);
    k_s(:,1,:) = k_s(:,1,:)/2;
    k_pad = ifftshift(padarray(padarray(k_s,[0 round(pad/2)],'pre'), [0 pad-round(pad/2)], 'post'),2);
end
img = ifft(ifft(k_pad,[],1),[],2);


% scanner frame
img = permute(img, [2 1 3 4]);
img = flipdim(flipdim(img,2),3);
[nv,np,ns,~] = size(img); % phase, readout, slice, receivers
vox = [Pars.lpe/nv, Pars.lro/np, Pars.lpe2/ns]*10;

% field directions
%% intrinsic euler angles 
% z-x-z convention, psi first, then theta, lastly phi
% psi and theta are left-handed, while gamma is right-handed!
alpha = - Pars.psi/180*pi;
beta = - Pars.theta/180*pi;
gamma =  Pars.phi/180*pi;
z_prjs = [sin(beta)*sin(gamma), sin(beta)*cos(gamma), cos(beta)]
if ~ isequal(z_prjs,[0 0 1])
    disp('This is angled slicing');
    pwd
end

% combine receivers
if Pars.RCVRS_ == 4
    % combine RF coils
    disp('--> combine RF rcvrs ...');
    img_cmb = coils_cmb(img,vox,ref_coi,eig_rad);
else  % single channel  
    img_cmb = img;
end

% save nifti
mkdir('combine');
nii = make_nii(abs(img_cmb),vox);
save_nii(nii,'combine/mag_cmb.nii');


% %% center k-space correction (readout direction)
% k = ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(img_cmb,1),2),[],1),[],2),1),2);
% [~,Ind] = max(abs(k(:)));
% Ix = ceil(mod(Ind,np*nv)/nv);

% % Apply phase ramp
% pix = np/2-Ix; % voxel shift
% ph_ramp = exp(-sqrt(-1)*2*pi*pix*(-1/2:1/np:1/2-1/np));
% img_cmb = img_cmb.* repmat(ph_ramp,[nv 1 ns]);

% save nifti
nii = make_nii(angle(img_cmb),vox);
save_nii(nii,'combine/ph_cmb.nii');

clear img;


%%% generate brain mask
disp('--> extract brain volume and generate mask ...');
setenv('bet_thr',num2str(bet_thr));
unix('bet combine/mag_cmb.nii BET -f ${bet_thr} -m -R');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);

if options.ero
    % erode the brain 2mm 
    imsize = size(mask);
    % make spherical/ellipsoidal convolution kernel (ker)
    rx = round(2/vox(1));
    ry = round(2/vox(2));
    rz = round(2/vox(3));
    % rz = ceil(ker_rad/vox(3));
    [X,Y,Z] = ndgrid(-rx:rx,-ry:ry,-rz:rz);
    h = (X.^2/rx^2 + Y.^2/ry^2 + Z.^2/rz^2 < 1);
    ker = h/sum(h(:));
    % circularshift, linear conv to Fourier multiplication
    csh = [rx,ry,rz]; % circularshift
    % erode the mask by convolving with the kernel
    cvsize = imsize + [2*rx+1, 2*ry+1, 2*rz+1] -1; % linear conv size
    mask_tmp = real(ifftn(fftn(mask,cvsize).*fftn(ker,cvsize)));
    mask_tmp = mask_tmp(rx+1:end-rx, ry+1:end-ry, rz+1:end-rz); % same size
    mask_ero = zeros(imsize);
    mask_ero(mask_tmp > 1-1/sum(h(:))) = 1;
    mask = mask_ero;
    unix('mv BET_mask.nii BET_mask_backup.nii');
    nii = make_nii(mask,vox);
    save_nii(nii,'BET_mask.nii');
end

% %% unwrap combined phase with PRELUDE
% disp('--> unwrap aliasing phase ...');
% unix('prelude -a combine/mag_cmb.nii -p combine/ph_cmb.nii -u unph.nii -m BET_mask.nii -n 8');
% unix('gunzip -f unph.nii.gz');
% nii = load_nii('unph.nii');
% unph = double(nii.img);


% % unwrap with Laplacian based method (TianLiu's)
% unph = unwrapLaplacian(angle(img_cmb), size(img_cmb), vox);
% nii = make_nii(unph, vox);
% save_nii(nii,'unph_lap.nii');


% Ryan Topfer's Laplacian unwrapping
Options.voxelSize = vox;
unph = lapunwrap(angle(img_cmb), Options);
nii = make_nii(unph, vox);
save_nii(nii,'unph_lap.nii');



% normalize to echo time and field strength
% ph = gamma*dB*TE
% dB/B = ph/(gamma*TE*B0)
% units: TE s, gamma 2.675e8 rad/(sT), B0 4.7T
% tfs = -unph_poly/(2.675e8*Pars.te*4.7)*1e6; % unit ppm
tfs = -unph/(2.675e8*Pars.te*4.7)*1e6; % unit ppm
nii = make_nii(tfs,vox);
save_nii(nii,'tfs.nii');




% (1) resharp 
disp('--> resharp to remove background field ...');
mkdir('resharp');
[lfs_resharp,mask_resharp] = resharp(tfs,mask,vox,smv_rad,tik_reg);

nii = make_nii(lfs_resharp,vox);
save_nii(nii,'resharp/lfs_resharp.nii');


%%% susceptibility inversion
disp('--> TV susceptibility inversion ...');
sus_resharp = tvdi(lfs_resharp,mask_resharp,vox,tv_reg,abs(img_cmb),z_prjs,tvdi_n);

nii = make_nii(sus_resharp,vox);
save_nii(nii,'resharp/sus_resharp.nii');

nii = make_nii(sus_resharp.*mask_resharp,vox);
save_nii(nii,'resharp/sus_resharp_clean.nii');



% (2) lbv
disp('--> lbv to remove background field ...');
lfs_lbv = LBV(tfs,mask,size(tfs),vox,0.01,2); % strip 2 layers
mkdir('lbv');
nii = make_nii(lfs_lbv,vox);
save_nii(nii,'lbv/lfs_lbv.nii');
mask_lbv = ones(size(mask));
mask_lbv(lfs_lbv==0) = 0;

% 3D 2nd order polyfit to remove phase-offset
lfs_lbv_poly= poly3d(lfs_lbv,mask_lbv);
nii = make_nii(lfs_lbv_poly,vox);
save_nii(nii,'lbv/lfs_lbv_poly.nii');


%%% susceptibility inversion
disp('--> TV susceptibility inversion ...');
% sus_lbv = tvdi(lfs_lbv,mask_lbv,vox,tv_reg,abs(img_cmb),z_prjs,tvdi_n);
sus_lbv = tvdi(lfs_lbv_poly,mask_lbv,vox,tv_reg,abs(img_cmb),z_prjs,tvdi_n);

nii = make_nii(sus_lbv,vox);
save_nii(nii,'lbv/sus_lbv.nii');

nii = make_nii(sus_lbv.*mask_lbv,vox);
save_nii(nii,'lbv/sus_lbv_clean.nii');



%%% save all variables for debugging purpose
if sav_all
    clear nii;
    save('all.mat','-v7.3');
end

% save parameters used in the recon
save('parameters.mat','options','-v7.3')


%%% clean up
% unix('rm *.nii*');
cd(init_dir);

