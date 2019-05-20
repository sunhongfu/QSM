function qsm_swi47(path_in, path_out, options)
%QSM_SWI47 Quantitative susceptibility mapping from SWI sequence at 4.7T.
%   QSM_SWI47(PATH_IN, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH_IN      - directory of .fid from ge3d sequence      : ge3d__01.fid
%   PATH_OUT     - directory to save nifti and/or matrixes   : QSM_SWI
%   OPTIONS      - parameter structure including fields below
%    .ref_coil   - reference coil to use for phase combine   : 2
%    .eig_rad    - radius (mm) of eig decomp kernel          : 5
%    .bet_thr    - threshold for BET brain mask              : 0.3
%    .bet_smooth - smoothness of BET brain mask at edges     : 2
%    .ph_unwrap  - 'prelude' or 'laplacian' or 'bestpath'    : 'laplacian'
%    .bkg_rm     - background field removal method(s)        : 'resharp'
%	               options: 'pdf','sharp','resharp','esharp','lbv'
%    .smv_rad    - radius (mm) of SMV convolution kernel     : 4
%    .tik_reg    - Tikhonov regularization for resharp       : 5e-4
%    .lbv_layer  - LBV layers to be stripped off             : 2
%    .t_svd      - truncation of SVD for SHARP               : 0.05
%    .tv_reg     - Total variation regularization parameter  : 5e-4
%    .tvdi_n     - iteration number of TVDI (nlcg)           : 500
%    .clean_all  - clean all the temp nifti results          : 1
%    .interp     - interpolate the image to the double size  : 0


% default settings
if ~ exist('path_in','var') || isempty(path_in)
    path_in = pwd;
end

if exist([path_in '/fid'],'file')
    path_fid = path_in;
    path_fid = cd(cd(path_fid));
elseif exist([path_in '/ge3d__01.fid/fid'],'file')
    path_fid = [path_in, '/ge3d__01.fid'];
    path_fid = cd(cd(path_fid));
else
    error('cannot find .fid file');
end

if ~ exist('path_out','var') || isempty(path_out)
    path_out = path_fid;
end

if ~ exist('options','var') || isempty(options)
    options = [];
end

if ~ isfield(options,'ref_coil')
    options.ref_coil = 2;
end

if ~ isfield(options,'eig_rad')
    options.eig_rad = 5;
end

if ~ isfield(options,'bet_thr')
    options.bet_thr = 0.3;
end

if ~ isfield(options,'bet_smooth')
    options.bet_smooth = 2;
end

if ~ isfield(options,'ph_unwrap')
    options.ph_unwrap = 'laplacian';
    % lap can handle singularities
end

if ~ isfield(options,'bkg_rm')
    options.bkg_rm = 'resharp';
    % options.bkg_rm = {'pdf','sharp','resharp','lbv'};
end

if ~ isfield(options,'t_svd')
    options.t_svd = 0.05;
end

if ~ isfield(options,'smv_rad')
    options.smv_rad = 4;
end

if ~ isfield(options,'tik_reg')
    options.tik_reg = 5e-4;
end

if ~ isfield(options,'lbv_layer')
    options.lbv_layer = 2;
end

if ~ isfield(options,'tv_reg')
    options.tv_reg = 5e-4;
end

if ~ isfield(options,'inv_num')
    options.inv_num = 500;
end

if ~ isfield(options,'clean_all')
    options.clean_all = 1;
end

if ~ isfield(options,'interp')
    options.interp = 0;
end

ref_coil   = options.ref_coil;
eig_rad    = options.eig_rad;
bet_thr    = options.bet_thr;
bet_smooth = options.bet_smooth;
ph_unwrap  = options.ph_unwrap;
bkg_rm     = options.bkg_rm;
t_svd      = options.t_svd;
smv_rad    = options.smv_rad;
tik_reg    = options.tik_reg;
lbv_layer  = options.lbv_layer;
tv_reg     = options.tv_reg;
inv_num    = options.inv_num;
clean_all  = options.clean_all;
interp     = options.interp;


%%% define directories
% if strcmpi(ph_unwrap,'prelude')
%     path_qsm = [path_out, filesep, 'QSM_SWI47_pre'];
% elseif strcmpi(ph_unwrap,'laplacian')
%     path_qsm = [path_out, filesep, 'QSM_SWI47_lap'];
% elseif strcmpi(ph_unwrap,'bestpath')
%     path_qsm = [path_out, filesep, 'QSM_SWI47_best'];
% end
path_qsm = [path_out, filesep, 'QSM_SWI47'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);


%%% generate raw img
disp('--> reconstruct fid to complex img ...');
% [img,Pars] = swi47_recon(path_fid,swi_ver);
[img,Pars] = swi47_recon(path_fid);


% match scanner frame (PE,RO,SL,NE,RX)
% so that angle corrections can be performed (phi, psi, theta)
img = permute(img, [2 1 3 4]);
img = flipdim(flipdim(img,2),3);
[nv,np,ns,nrcvrs] = size(img); % phase, readout, slice, receivers
voxelSize = [Pars.lpe/nv, Pars.lro/np, Pars.lpe2/ns]*10;


% have a peak of the raw phase
nii = make_nii(angle(img),voxelSize);
save_nii(nii,'rawphase.nii');


% intrinsic euler angles 
% z-x-z convention, psi first, then theta, lastly phi
% psi and theta are left-handed, while gamma is right-handed!
% alpha = - Pars.psi/180*pi;
beta = - Pars.theta/180*pi;
gamma =  Pars.phi/180*pi;
z_prjs = [sin(beta)*sin(gamma), sin(beta)*cos(gamma), cos(beta)];
if ~ isequal(z_prjs,[0 0 1])
    disp('This is angled slicing');
    disp(z_prjs);
    pwd
end


% center k-space correction (readout direction)
ks = ifftshift(ifftn(img));
[MAX,Ind] = max(abs(ks(:)));;
% find maximum readout and phase encoding index
[Inv, Inp, Ins, Ircvrs] = ind2sub([nv,np,ns,nrcvrs],Ind);

% Apply phase ramp
pix = np/2-Inp; % voxel shift
ph_ramp = exp(-sqrt(-1)*2*pi*pix*(-1/2:1/np:1/2-1/np));
pix2 = nv/2-Inv;
ph_ramp2 = exp(-sqrt(-1)*2*pi*pix2*(-1/2:1/nv:1/2-1/nv));

img_corr = img.* repmat((ph_ramp),[nv 1 ns nrcvrs]);
img_corr = img_corr.* repmat(transpose(ph_ramp2),[1 np ns nrcvrs]);


% combine receivers
if Pars.RCVRS_ > 1
    % combine RF coils
    disp('--> combine RF rcvrs ...');
    img_cmb = adaptive_cmb(img_corr,voxelSize,ref_coil,eig_rad);
else  % single channel  
    img_cmb = img_corr;
end

% interpolate the images to the double size
if interp
    imsize = size(img_cmb);
    img_cmb = single(img_cmb);
    % zero padding the k-space
    k = fftshift(fftshift(fftshift(fft(fft(fft(img_cmb,[],1),[],2),[],3),1),2),3);
    k = padarray(k,double([0 0 25]));
    img_cmb = ifft(ifft(ifft(ifftshift(ifftshift(ifftshift(k,1),2),3),[],1),[],2),[],3);
    clear k;
    imsize = size(img_cmb);
    voxelSize(3) = voxelSize(3)/2;
end


% save nifti
mkdir('combine');
nii = make_nii(abs(img_cmb),voxelSize);
save_nii(nii,'combine/mag_cmb.nii');
nii = make_nii(angle(img_cmb),voxelSize);
save_nii(nii,'combine/ph_cmb.nii');


% generate brain mask
disp('--> extract brain volume and generate mask ...');
setenv('bet_thr',num2str(bet_thr));
setenv('bet_smooth',num2str(bet_smooth));
[status,cmdout] = unix('rm BET*');
unix('bet2 combine/mag_cmb.nii BET -f ${bet_thr} -m -w ${bet_smooth}');
% unix('bet2 combine/mag_cmb.nii BET -f ${bet_thr} -m');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);


% phase unwrapping, prelude is preferred!
if strcmpi('prelude',ph_unwrap)
    % unwrap combined phase with PRELUDE
    disp('--> unwrap aliasing phase using prelude...');
    bash_script = ['prelude -a combine/mag_cmb.nii -p combine/ph_cmb.nii ' ...
        '-u unph.nii -m BET_mask.nii -n 12'];
    unix(bash_script);
    unix('gunzip -f unph.nii.gz');
    nii = load_nii('unph.nii');
    unph = double(nii.img);

    % % unwrap with Laplacian based method (TianLiu's)
    % unph = unwrapLaplacian(angle(img_cmb), size(img_cmb), voxelSize);
    % nii = make_nii(unph, voxelSize);
    % save_nii(nii,'unph_lap.nii');

elseif strcmpi('laplacian',ph_unwrap)
    % Ryan Topfer's Laplacian unwrapping
    disp('--> unwrap aliasing phase using laplacian...');
    Options.voxelSize = voxelSize;
    unph = lapunwrap(angle(img_cmb), Options);
    nii = make_nii(unph, voxelSize);
    save_nii(nii,'unph_lap.nii');

elseif strcmpi('bestpath',ph_unwrap)
    % unwrap the phase using best path
    disp('--> unwrap aliasing phase using bestpath...');
        [pathstr, ~, ~] = fileparts(which('3DSRNCP.m'));
    setenv('pathstr',pathstr);
    setenv('nv',num2str(nv));
    setenv('np',num2str(np));
    setenv('ns',num2str(ns));

    fid = fopen('wrapped_phase.dat','w');
    fwrite(fid,angle(img_cmb),'float');
    fclose(fid);
    mask_unwrp = uint8(abs(mask)*255);
    fid = fopen('mask_unwrp.dat','w');
    fwrite(fid,mask_unwrp,'uchar');
    fclose(fid);

    bash_script = ['${pathstr}/3DSRNCP wrapped_phase.dat mask_unwrp.dat unwrapped_phase.dat ' ...
        '$nv $np $ns reliability.dat'];
    unix(bash_script);

    fid = fopen('unwrapped_phase.dat','r');
    tmp = fread(fid,'float');
    unph = reshape(tmp - round(mean(tmp(mask==1))/(2*pi))*2*pi,[nv,np,ns]).*mask;
    fclose(fid);

    nii = make_nii(unph,voxelSize);
    save_nii(nii,'unph_bestpath.nii');

    fid = fopen('reliability.dat','r');
    reliability_raw = fread(fid,'float');
    reliability_raw = reshape(reliability_raw,[nv,np,ns]);
    fclose(fid);

    nii = make_nii(reliability_raw.*mask,voxelSize);
    save_nii(nii,'reliability_raw.nii');
    
    % reliability = mask;
    % reliability(reliability_raw >= 20) = 0;
    % % reliability(reliability > 0.1) = 1;
    % nii = make_nii(reliability,voxelSize);
    % save_nii(nii,'reliability.nii');
    % weights = abs(img_cmb)./max(abs(img_cmb(:))).*mask.*reliability;
    % weights = smooth3(weights,'gaussian',[7,7,3],1);
    % % weights = smooth3(weights,'gaussian',[7,7,3],0.5);
    % nii = make_nii(weights,voxelSize);
    % save_nii(nii,'weights.nii');
else
    error('what unwrapping methods to use? prelude or laplacian or bestpath?')
end


% normalize to echo time and field strength
% ph = gamma*dB*TE
% dB/B = ph/(gamma*TE*B0)
% units: TE s, gamma 2.675e8 rad/(sT), B0 4.7T
% tfs = -unph_poly/(2.675e8*Pars.te*4.7)*1e6; % unit ppm
tfs = -unph/(2.675e8*Pars.te*4.7)*1e6; % unit ppm
nii = make_nii(tfs,voxelSize);
save_nii(nii,'tfs.nii');


if ~ exist('weights')
    weights = abs(img_cmb);
end

% PDF
if sum(strcmpi('pdf',bkg_rm))
    disp('--> PDF to remove background field ...');
    [lfs_pdf,mask_pdf] = pdf(tfs,mask,voxelSize,smv_rad, ...
        weights,z_prjs);
    % 3D 2nd order polyfit to remove any residual background
    lfs_pdf= lfs_pdf - poly3d(lfs_pdf,mask_pdf);

    % save nifti
    mkdir('PDF');
    nii = make_nii(lfs_pdf,voxelSize);
    save_nii(nii,'PDF/lfs_pdf.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on PDF...');
    sus_pdf = tvdi(lfs_pdf, mask_pdf, voxelSize, tv_reg, ...
        weights, z_prjs, inv_num); 

    % save nifti
    nii = make_nii(sus_pdf.*mask_pdf,voxelSize);
    save_nii(nii,'PDF/sus_pdf.nii');
end


% SHARP (t_svd: truncation threthold for t_svd)
if sum(strcmpi('sharp',bkg_rm))
    disp('--> SHARP to remove background field ...');
    [lfs_sharp, mask_sharp] = sharp(tfs,mask,voxelSize,smv_rad,t_svd);
    % 3D 2nd order polyfit to remove any residual background
    lfs_sharp= lfs_sharp - poly3d(lfs_sharp,mask_sharp);

    % save nifti
    mkdir('SHARP');
    nii = make_nii(lfs_sharp,voxelSize);
    save_nii(nii,'SHARP/lfs_sharp.nii');
    
    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on SHARP...');
    sus_sharp = tvdi(lfs_sharp, mask_sharp, voxelSize, tv_reg, ...
        weights, z_prjs, inv_num); 
   
    % save nifti
    nii = make_nii(sus_sharp.*mask_sharp,voxelSize);
    save_nii(nii,'SHARP/sus_sharp.nii');
end


% RE-SHARP (tik_reg: Tikhonov regularization parameter)
if sum(strcmpi('resharp',bkg_rm))
    disp('--> RESHARP to remove background field ...');
    [lfs_resharp, mask_resharp] = resharp(tfs,mask,voxelSize,smv_rad,tik_reg);
    % 3D 2nd order polyfit to remove any residual background
    lfs_resharp= lfs_resharp - poly3d(lfs_resharp,mask_resharp);

    % save nifti
    mkdir('RESHARP');
    nii = make_nii(lfs_resharp,voxelSize);
    save_nii(nii,'RESHARP/lfs_resharp.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on RESHARP...');
    sus_resharp = tvdi(lfs_resharp, mask_resharp, voxelSize, tv_reg, ...
        weights, z_prjs, inv_num); 
   
    % save nifti
    nii = make_nii(sus_resharp.*mask_resharp,voxelSize);
    save_nii(nii,'RESHARP/sus_resharp.nii');
end


% E-SHARP (SHARP edge extension)
if sum(strcmpi('esharp',bkg_rm))
    disp('--> E-SHARP to remove background field ...');
    Parameters.voxelSize             = voxelSize; % in mm
    Parameters.resharpRegularization = tik_reg ;
    Parameters.resharpKernelRadius   = smv_rad ; % in mm
    Parameters.radius                = [ 10 10 5 ] ;

% pad matrix size to even number
    pad_size = mod(size(tfs),2);
    tfs = padarray(tfs, pad_size, 'post');
    mask = padarray(mask, pad_size, 'post');

    % taking off additional 3 voxels from edge - not sure the outermost 
    % phase data included in the original mask is reliable. 
    tfs        = tfs .* mask;
    mask       = shaver( ( tfs ~= 0 ), 1 ) ; % 1 voxel taken off
    totalField = mask .* tfs ;

    % resharp 
    [reducedLocalField, maskReduced] = ...
        resharp( totalField, ...
                 double(mask), ...
                 Parameters.voxelSize, ...
                 Parameters.resharpKernelRadius, ...
                 Parameters.resharpRegularization ) ;

    % extrapolation ~ esharp 
    reducedBackgroundField = maskReduced .* ( totalField - reducedLocalField) ;

    extendedBackgroundField = extendharmonicfield( ...
       reducedBackgroundField, mask, maskReduced, Parameters) ;

    backgroundField = extendedBackgroundField + reducedBackgroundField ;
    localField      = totalField - backgroundField ;

    lfs_esharp      = localField(1+pad_size(1):end,1+pad_size(2):end,1+pad_size(3):end);
    mask_esharp     = mask(1+pad_size(1):end,1+pad_size(2):end,1+pad_size(3):end);  

    % 3D 2nd order polyfit to remove any residual background
    lfs_esharp = lfs_esharp - poly3d(lfs_esharp,mask_esharp);

    % save nifti
    mkdir('ESHARP');
    nii = make_nii(lfs_esharp,voxelSize);
    save_nii(nii,'ESHARP/lfs_esharp.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on ESHARP...');
    sus_esharp = tvdi(lfs_esharp, mask_esharp, voxelSize, tv_reg, ...
        weights, z_prjs, inv_num); 
   
    % save nifti
    nii = make_nii(sus_esharp.*mask_esharp,voxelSize);
    save_nii(nii,'ESHARP/sus_esharp.nii');
end


% LBV
if sum(strcmpi('lbv',bkg_rm))
    disp('--> LBV to remove background field ...');
    lfs_lbv = LBV(tfs,mask,size(tfs),voxelSize,0.01,lbv_layer); % strip 2 layers
    mask_lbv = ones(size(mask));
    mask_lbv(lfs_lbv==0) = 0;
    % 3D 2nd order polyfit to remove any residual background
    lfs_lbv= lfs_lbv - poly3d(lfs_lbv,mask_lbv);

    % save nifti
    mkdir('LBV');
    nii = make_nii(lfs_lbv.*mask_lbv,voxelSize);
    save_nii(nii,'LBV/lfs_lbv.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on lbv...');
    sus_lbv = tvdi(lfs_lbv,mask_lbv,voxelSize,tv_reg, ...
        weights,z_prjs,inv_num);   

    % save nifti
    nii = make_nii(sus_lbv.*mask_lbv,voxelSize);
    save_nii(nii,'LBV/sus_lbv.nii');
end


% clean the directory
if clean_all
    disp('--> clean temp nifti files ...');
    unix('ls | grep -v "combine\|RESHARP\|unph" | xargs rm -rf');
else
    % save all variables for future reference
    clear nii;
    disp('--> save the entire workspace ...');
    save('all.mat','-v7.3');
end

% save parameters used in the recon
save('parameters.mat','options','-v7.3')

% save the git log for future tracking
unix('git log --branches --decorate --color --abbrev-commit --graph --no-merges --tags > git_log');

% save all the NIFTIs in LPS orientation
% originally in scanner coordinates LAI
mkdir('LPS');
[status, list] = unix('find . -name "*.nii"');
expression = '\n';
splitStr = regexp(strtrim(list),expression,'split');
for i = 1:size(splitStr,2)
    [pathstr,name,ext] = fileparts(splitStr{i});
    nii = load_nii(splitStr{i});
    tmp = double(nii.img);
    tmp = flipdim(flipdim(tmp,2),3);
    nii = make_nii(tmp,voxelSize);
    save_nii(nii,['LPS/',name,ext]);
end

% go back to the initial directory
cd(init_dir);
