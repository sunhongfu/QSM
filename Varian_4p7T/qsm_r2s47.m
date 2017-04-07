function qsm_r2s47(path_in, path_out, options)
%QSM_R2S47 Quantitative susceptibility mapping from R2* sequence at 4.7T.
%   QSM_R2S47(PATH_IN, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH_IN      - directory of .fid from gemsme3d sequence  : gemsme3d_R2s_01.fid
%   PATH_OUT     - directory to save nifti and/or matrixes   : QSM_R2s47
%   OPTIONS      - parameter structure including fields below
%    .ref_coil   - reference coil to use for phase combine   : 2
%    .eig_rad    - radius (mm) of eig decomp kernel          : 5
%    .bet_thr    - threshold for BET brain mask              : 0.5
%    .bet_smooth - smoothness of BET brain mask at edges     : 2
%    .bkg_rm     - background field removal method(s)        : 'resharp'
%	               choices: 'pdf','sharp','resharp','esharp','lbv'
%    .smv_rad    - radius (mm) of SMV convolution kernel     : 4
%    .tik_reg    - Tikhonov regularization for RESHARP       : 0.0005
%    .t_svd      - truncation of SVD for SHARP               : 0.05
%    .tv_reg     - Total variation regularization parameter  : 0.0005
%    .inv_num    - iteration number of TVDI (nlcg)           : 500
%    .echo_num   - keep only the first 'echo_num' echoes     : 5
%    .fit_thr    - truncation level on fitting residual      : 10
%    .clean_all  - clean all the temp nifti results          : 1
%    .interp     - interpolate the image to the double size  : 0


% default settings
if ~ exist('path_in','var') || isempty(path_in)
    path_in = pwd;
end

if exist([path_in '/fid'],'file')
    path_fid = path_in;
    path_fid = cd(cd(path_fid));
elseif exist([path_in '/gemsme3d_R2s_01.fid/fid'],'file')
    path_fid = [path_in, '/gemsme3d_R2s_01.fid'];
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
    options.bet_thr = 0.5;
end

if ~ isfield(options,'bet_smooth')
    options.bet_smooth = 2;
end

if ~ isfield(options,'echo_num')
    options.echo_num = 5;
end

if ~ isfield(options,'r_mask')
    options.r_mask = 1;
end

if ~ isfield(options,'fit_thr')
    options.fit_thr = 10;
end

if ~ isfield(options,'bkg_rm')
    options.bkg_rm = 'resharp';
    % options.bkg_rm = {'pdf','sharp','resharp','lbv'};
end

if ~ isfield(options,'smv_rad')
    options.smv_rad = 4;
end

if ~ isfield(options,'tik_reg')
    options.tik_reg = 5e-4;
end

if ~ isfield(options,'t_svd')
    options.t_svd = 0.05;
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
echo_num   = options.echo_num;
r_mask     = options.r_mask; 
fit_thr    = options.fit_thr;
bkg_rm     = options.bkg_rm;
smv_rad    = options.smv_rad;
tik_reg    = options.tik_reg;
t_svd      = options.t_svd;
tv_reg     = options.tv_reg;
inv_num    = options.inv_num;
clean_all  = options.clean_all;
interp     = options.interp;


% define directories
path_qsm = [path_out '/QSM_R2s47'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);


% reconstruct complex image from fid file
disp('--> reconstruct fid to complex img ...');
[img,par] = r2s47_recon(path_fid);

imsize = size(img);

% interpolate the images to the double size
if interp
    img = single(img);
    % zero padding the k-space
    k = fftshift(fftshift(fftshift(fft(fft(fft(img,[],1),[],2),[],3),1),2),3);
    k = padarray(k,double(imsize(1:3)/2));
    img = ifft(ifft(ifft(ifftshift(ifftshift(ifftshift(k,1),2),3),[],1),[],2),[],3);
    clear k;
    imsize = size(img);
    % vox = vox/2;
end


% match scanner frame (PE,RO,SL,NE,RX)
% so that angle corrections can be performed (phi, psi, theta)
img = permute(img,[2 1 3 4 5]); 
img = flipdim(img,1);
img = flipdim(img,2);
[nv,np,nv2,ne,~] = size(img);
voxelSize = [par.lpe/nv, par.lro/np, par.lpe2/nv2]*10;
% resolution in mm/pixel
te = par.te + (0:ne-1)*par.esp;


% intrinsic euler angles 
% z-x-z convention, psi first, then theta, lastly phi
% psi and theta are left-handed, while gamma is right-handed!
% alpha = - par.psi/180*pi;
beta = - par.theta/180*pi;
gamma =  par.phi/180*pi;
z_prjs = [sin(beta)*sin(gamma), sin(beta)*cos(gamma), cos(beta)];
if ~ isequal(z_prjs,[0 0 1])
    disp('This is not pure axial slicing');
    disp(z_prjs);
end

% have a peak at raw phase
nii = make_nii(squeeze(angle(img(:,:,:,1,:))));
save_nii(nii,'rawphase.nii');

% combine magnitudes using eig method (DO Walsh, MRM2000)
if par.nrcvrs > 1
    disp('--> combine RF rcvrs ...');
    img_cmb = adaptive_cmb(permute(img,[1 2 3 5 4]),voxelSize,ref_coil,eig_rad);
    mag_cmb = abs(img_cmb);
    % at 4.7T, seems the 2rd coil has the best SNR?
else
    mag_cmb = abs(img);
end

% save niftis after coil combination
mkdir('combine');
for echo = 1:ne
    nii = make_nii(mag_cmb(:,:,:,echo),voxelSize);
    save_nii(nii,['combine/mag_cmb' num2str(echo) '.nii']);
end


% generate mask from combined magnitude of the 1th echo
disp('--> extract brain volume and generate mask ...');
setenv('bet_thr',num2str(bet_thr));
setenv('bet_smooth',num2str(bet_smooth));
[status,cmdout] = unix('rm BET*');
% unix('bet2 combine/mag_cmb1.nii BET -f ${bet_thr} -m -w ${bet_smooth}');
unix('bet2 combine/mag_cmb1.nii BET -f ${bet_thr} -m');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);



% combine phase using double-echo method
% always use geme_cmb even for only 1 receiver
% this function can properly remove offset
% if par.nrcvrs > 1
    ph_cmb = geme_cmb(img,voxelSize,te,mask);
% else
%     ph_cmb = angle(img);
% end


% save niftis after coil combination
for echo = 1:size(ph_cmb,4)
    nii = make_nii(ph_cmb(:,:,:,echo),voxelSize);
    save_nii(nii,['combine/ph_cmb' num2str(echo) '.nii']);
end

img_cmb = mag_cmb.*exp(1j.*ph_cmb);
% clear mag_cmb ph_cmb


% keep only the first 'echo_num' echoes
echo_num = min(ne,echo_num);
img_cmb = img_cmb(:,:,:,1:echo_num);
ne = echo_num;
te = par.te + (0:ne-1)*par.esp;


% unwrap phase from each echo
disp('--> unwrap aliasing phase for all TEs using prelude...');

setenv('echo_num',num2str(echo_num));
bash_command = sprintf(['for ph in combine/ph_cmb[1-$echo_num].nii\n' ...
'do\n' ...
'	base=`basename $ph`;\n' ...
'	dir=`dirname $ph`;\n' ...
'	mag=$dir/"mag"${base:2};\n' ...
'	unph="unph"${base:2};\n' ...
'	prelude -a $mag -p $ph -u $unph -m BET_mask.nii -n 12&\n' ...
'done\n' ...
'wait\n' ...
'gunzip -f unph*.gz\n']);

unix(bash_command);

unph_cmb = zeros(nv,np,nv2,ne);
for echo = 1:ne
    nii = load_nii(['unph_cmb' num2str(echo) '.nii']);
    unph_cmb(:,:,:,echo) = double(nii.img);
end


% check and correct for 2pi jump between echoes
disp('--> correct for potential 2pi jumps between TEs ...')

% nii = load_nii('unph_cmb1.nii');
% unph1 = double(nii.img);
% nii = load_nii('unph_cmb2.nii');
% unph2 = double(nii.img);
% unph_diff = unph2 - unph1;

nii = load_nii('unph_diff.nii');
unph_diff = double(nii.img);

for echo = 2:ne
    meandiff = unph_cmb(:,:,:,echo)-unph_cmb(:,:,:,1)-(echo-1)*unph_diff;
    meandiff = meandiff(mask==1);
    meandiff = mean(meandiff(:))
    njump = round(meandiff/(2*pi))
    disp(['    ' num2str(njump) ' 2pi jumps for TE' num2str(echo)]);
    unph_cmb(:,:,:,echo) = unph_cmb(:,:,:,echo) - njump*2*pi;
    unph_cmb(:,:,:,echo) = unph_cmb(:,:,:,echo).*mask;
end


% fit phase images with echo times
disp('--> magnitude weighted LS fit of phase to TE ...');
[tfs, fit_residual] = echofit(unph_cmb,abs(img_cmb),te); 

% normalize to main field
% ph = gamma*dB*TE
% dB/B = ph/(gamma*TE*B0)
% units: TE s, gamma 2.675e8 rad/(sT), B0 4.7T
tfs = -tfs/(2.675e8*4.7)*1e6; % unit ppm


if r_mask
    % generate reliability map
    fit_residual_blur = smooth3(fit_residual,'box',round(smv_rad./voxelSize/2)*2+1); 
    nii = make_nii(fit_residual_blur,voxelSize);
    save_nii(nii,'fit_residual_blur.nii');
    R = ones(size(fit_residual_blur));
    R(fit_residual_blur >= fit_thr) = 0;
else
    R = 1;
end


% PDF
if sum(strcmpi('pdf',bkg_rm))
    disp('--> PDF to remove background field ...');
    [lfs_pdf,mask_pdf] = pdf(tfs,mask.*R,voxelSize,smv_rad, ...
        abs(img_cmb(:,:,:,end)),z_prjs);
    % % 3D 2nd order polyfit to remove any residual background
    % lfs_pdf= poly3d(lfs_pdf,mask_pdf);

    % save nifti
    mkdir('PDF');
    nii = make_nii(lfs_pdf,voxelSize);
    save_nii(nii,'PDF/lfs_pdf.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on PDF...');
    sus_pdf = tvdi(lfs_pdf, mask_pdf, voxelSize, tv_reg, ...
        abs(img_cmb(:,:,:,echo_num)), z_prjs, inv_num); 

    % save nifti
    nii = make_nii(sus_pdf.*mask_pdf,voxelSize);
    save_nii(nii,'PDF/sus_pdf.nii');
end


% SHARP (t_svd: truncation threthold for t_svd)
if sum(strcmpi('sharp',bkg_rm))
    disp('--> SHARP to remove background field ...');
    [lfs_sharp, mask_sharp] = sharp(tfs,mask.*R,voxelSize,smv_rad,t_svd);
    % % 3D 2nd order polyfit to remove any residual background
    % lfs_sharp= poly3d(lfs_sharp,mask_sharp);

    % save nifti
    mkdir('SHARP');
    nii = make_nii(lfs_sharp,voxelSize);
    save_nii(nii,'SHARP/lfs_sharp.nii');
    
    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on SHARP...');
    sus_sharp = tvdi(lfs_sharp, mask_sharp, voxelSize, tv_reg, ...
        abs(img_cmb(:,:,:,echo_num)), z_prjs, inv_num); 
   
    % save nifti
    nii = make_nii(sus_sharp.*mask_sharp,voxelSize);
    save_nii(nii,'SHARP/sus_sharp.nii');
end


% RE-SHARP (tik_reg: Tikhonov regularization parameter)
if sum(strcmpi('resharp',bkg_rm))
    disp('--> RESHARP to remove background field ...');
    [lfs_resharp, mask_resharp] = resharp(tfs,mask.*R,voxelSize,smv_rad,tik_reg);


    % save nifti
    mkdir('RESHARP');
    nii = make_nii(lfs_resharp,voxelSize);
    save_nii(nii,'RESHARP/lfs_resharp.nii');


    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on RESHARP...');
    sus_resharp = tvdi(lfs_resharp, mask_resharp, voxelSize, tv_reg, ...
        abs(img_cmb(:,:,:,echo_num)), z_prjs, inv_num); 
   

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
    tfs = double(padarray(tfs, pad_size, 'post'));
    mask = double(padarray(mask, pad_size, 'post'));

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

    % save nifti
    mkdir('ESHARP');
    nii = make_nii(lfs_esharp,voxelSize);
    save_nii(nii,'ESHARP/lfs_esharp.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on ESHARP...');
    sus_esharp = tvdi(lfs_esharp, mask_esharp, voxelSize, tv_reg, ...
        abs(img_cmb(:,:,:,echo_num)), z_prjs, inv_num);

    % save nifti
    nii = make_nii(sus_esharp.*mask_esharp,voxelSize);
    save_nii(nii,'ESHARP/sus_esharp.nii');
end




% LBV
if sum(strcmpi('lbv',bkg_rm))
    disp('--> LBV to remove background field ...');
    lfs_lbv = LBV(tfs,mask.*R,size(tfs),voxelSize,0.01,2); % strip 2 layers
    mask_lbv = ones(size(mask));
    mask_lbv(lfs_lbv==0) = 0;
    % % 3D 2nd order polyfit to remove phase-offset
    % lfs_lbv= poly3d(lfs_lbv,mask_lbv);


    % save nifti
    mkdir('LBV');
    nii = make_nii(lfs_lbv,voxelSize);
    save_nii(nii,'LBV/lfs_lbv.nii');


    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on lbv...');
    sus_lbv = tvdi(lfs_lbv,mask_lbv,voxelSize,tv_reg, ...
        abs(img_cmb(:,:,:,echo_num)),z_prjs,inv_num);   

    % save nifti
    nii = make_nii(sus_lbv.*mask_lbv,voxelSize);
    save_nii(nii,'LBV/sus_lbv.nii');
end



% clean the directory
if clean_all
    disp('--> clean temp nifti files ...');
    unix('ls | grep -v "combine\|RESHARP" | xargs rm -rf');
else
    % save all variables for future reference
    clear nii;
    disp('--> save the entire workspace ...');
    save('all.mat','-v7.3');
end

% save parameters used in the recon
save('parameters.mat','options','-v7.3');

% save the git log for future tracking
unix('git log --branches --decorate --color --abbrev-commit --graph --no-merges --tags > git_log');

% go back to the initial directory
cd(init_dir);


