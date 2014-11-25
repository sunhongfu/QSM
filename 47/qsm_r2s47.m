function qsm_r2s47(path_in, path_out, options)
%QSM_R2S47 Quantitative susceptibility mapping from R2* sequence at 4.7T.
%   QSM_R2S47(PATH_IN, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH_IN    - directory of .fid from gemsme3d sequence  : gemsme3d_R2s_01.fid
%   PATH_OUT   - directory to save nifti and/or matrixes   : QSM_R2s_vxxx
%   OPTIONS    - parameter structure including fields below
%    .bet_thr  - threshold for BET brain mask              : 0.5
%    .bkg_rm    - background field removal method(s)        : 'resharp'
%    .smv_rad  - radius (mm) of SMV convolution kernel     : 4
%    .tik_reg  - Tikhonov regularization for RESHARP       : 0.0005
%    .t_svd     - truncation of SVD for SHARP               : 0.05
%    .tv_reg   - Total variation regularization parameter  : 0.0005
%    .inv_num   - iteration number of TVDI (nlcg)           : 200
%    .echo_num   - keep only the first 'echo_num' echoes       : 5
%    .fit_thr    - truncation level on fitting residual      : 10
%    .save_all  - save all the variables for debug (~ 0)    : 1


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
    options.ref_coil = 3;
end

if ~ isfield(options,'eig_rad')
    options.eig_rad = 4;
end

if ~ isfield(options,'bet_thr')
    options.bet_thr = 0.5;
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
    % options.bkg_rm = {'resharp','lbv'};
    % options.bkg_rm = 'resharp';
    options.bkg_rm = {'pdf','sharp','resharp','lbv'};
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
    options.inv_num = 200;
end

if ~ isfield(options,'save_all')
    options.save_all = 1;
end


ref_coil = options.ref_coil;
eig_rad  = options.eig_rad;
bet_thr  = options.bet_thr;
echo_num = options.echo_num;
r_mask   = options.r_mask; 
fit_thr  = options.fit_thr;
bkg_rm   = options.bkg_rm;
smv_rad  = options.smv_rad;
tik_reg  = options.tik_reg;
t_svd    = options.t_svd;
tv_reg   = options.tv_reg;
inv_num  = options.inv_num;
save_all = options.save_all;


% define directories
path_qsm = [path_out '/QSM_R2s_v500'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);




% reconstruct complex image from fid file
disp('--> reconstruct fid to complex img ...');
[img,par] = r2s47_recon(path_fid);


% match scanner frame (PE,RO,SL,NE,RX)
% so that angle corrections can be performed (phi, psi, theta)
img = permute(img,[2 1 3 4 5]); 
img = flipdim(img,1);
img = flipdim(img,2);
[nv,np,nv2,ne,~] = size(img)
voxelSize = [par.lpe/par.nv, par.lro/(par.np/2), par.lpe2/par.nv2]*10;
% resolution in mm/pixel
te = par.te + (0:ne-1)*par.esp;


% intrinsic euler angles 
% z-x-z convention, psi first, then theta, lastly phi
% psi and theta are left-handed, while gamma is right-handed!
alpha = - par.psi/180*pi;
beta = - par.theta/180*pi;
gamma =  par.phi/180*pi;
z_prjs = [sin(beta)*sin(gamma), sin(beta)*cos(gamma), cos(beta)];
if ~ isequal(z_prjs,[0 0 1])
    disp('This is not pure axial slicing');
    z_prjs
end


% combine magnitudes using eig method (DO Walsh, MRM2000)
if par.nrcvrs > 1
    img_cmb = coils_cmb(permute(img,[1 2 3 5 4]),voxelSize,ref_coil,eig_rad);
    mag_cmb = abs(img_cmb);
    % at 4.7T, seems the 3rd coil has the best SNR?
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
unix('bet combine/mag_cmb1.nii BET -f ${bet_thr} -m -R');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);


if ne > 1
    % combine phase using double-echo method
    if par.nrcvrs > 1
        ph_cmb = sense_me(img,voxelSize,te);
    else
        ph_cmb = angle(img);
    end
else
    if par.nrcvrs > 1
        ph_cmb = angle(img_cmb);
    else
        ph_cmb = angle(img);
    end
end

% save niftis after coil combination
mkdir('combine');
for echo = 1:size(ph_cmb,4)
    nii = make_nii(ph_cmb(:,:,:,echo),voxelSize);
    save_nii(nii,['combine/ph_cmb' num2str(echo) '.nii']);
end

img_cmb = mag_cmb.*exp(1j.*ph_cmb);
clear mag_cmb ph_cmb


% keep only the first 'echo_num' echoes
echo_num = min(ne,echo_num);
img_cmb = img_cmb(:,:,:,1:echo_num);
ne = echo_num;
te = par.te + (0:ne-1)*par.esp;




%% unwrap phase from each echo
disp('--> unwrap aliasing phase for all TEs ...');

setenv('echo_num',num2str(echo_num));
bash_command = sprintf(['for ph in combine/ph_cmb[1-$echo_num].nii\n' ...
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

unph_cmb = zeros(nv,np,nv2,ne);
for echo = 1:ne
    nii = load_nii(['unph_cmb' num2str(echo) '.nii']);
    unph_cmb(:,:,:,echo) = double(nii.img);
end


% check and correct for 2pi jump between echoes
disp('--> correct for potential 2pi jumps between TEs ...')

nii = load_nii('unph_cmb1.nii');
unph1 = double(nii.img);
nii = load_nii('unph_cmb2.nii');
unph2 = double(nii.img);
unph_diff = unph2 - unph1;

for echo = 2:ne
    meandiff = unph_cmb(:,:,:,echo)-unph_cmb(:,:,:,1)-(echo-1)*unph_diff;
    meandiff = meandiff(mask==1);
    meandiff = mean(meandiff(:));
    njump = round(meandiff/(2*pi));
    disp(['    ' num2str(njump) ' 2pi jumps for TE' num2str(echo)]);
    unph_cmb(:,:,:,echo) = unph_cmb(:,:,:,echo) - njump*2*pi;
    unph_cmb(:,:,:,echo) = unph_cmb(:,:,:,echo).*mask;
end


%% fit phase images with echo times
disp('--> magnitude weighted LS fit of phase to TE ...');
[tfs, fit_residual] = echofit(unph_cmb,abs(img_cmb),te); 

% normalize to main field
% ph = gamma*dB*TE
% dB/B = ph/(gamma*TE*B0)
% units: TE s, gamma 2.675e8 rad/(sT), B0 4.7T
tfs = -tfs/(2.675e8*4.7)*1e6; % unit ppm


if r_mask
    % generate reliability map
    fit_residual_blur = smooth3(fit_residual,'box',round(smv_rad./voxelSize)*2+1); 
    nii = make_nii(fit_residual_blur,voxelSize);
    save_nii(nii,'fit_residual_blur.nii');
    R = ones(size(fit_residual_blur));
    R(fit_residual_blur >= fit_thr) = 0;
else
    R = 1;
end


%% PDF
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


%% SHARP (t_svd: truncation threthold for t_svd)
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


%% RE-SHARP (tik_reg: Tikhonov regularization parameter)
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


%% LBV
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

%% E-SHARP
% if sum(strcmpi('esharp',bkg_rm))
%     disp('--> E-SHARP to remove background field ...');
%     Options.voxelSize = voxelSize;
%     lfs = esharp(tfs,mask.*R,Options);
%     mask_final = mask.*R;
% 
% 
%     % save nifti
%     mkdir([path_qsm '/rmbkg/']);
%     nii = make_nii(lfs,voxelSize);
%     save_nii(nii,[path_qsm '/rmbkg/lfs_esharp_xy.nii']);
%     nii = make_nii(permute(lfs,[1 3 2]),voxelSize);
%     save_nii(nii,[path_qsm '/rmbkg/lfs_esharp_xz.nii']);
%     nii = make_nii(permute(lfs,[2 3 1]),voxelSize);
%     save_nii(nii,[path_qsm '/rmbkg/lfs_esharp_yz.nii']);
%     nii = make_nii(mask_final,voxelSize);
%     save_nii(nii,[path_qsm '/mask/mask_esharp_final.nii']);
% 
% 
%     % inversion of susceptibility 
%     disp('--> (9/9) TV susceptibility inversion on E-SHARP...');
%     sus = tvdi(lfs, mask_final, voxelSize, tv_reg, mag_cmb(:,:,:,4)); 
% 
% 
%     % save nifti
%     mkdir([path_qsm '/inversion']);
%     nii = make_nii(sus,voxelSize);
%     save_nii(nii,[path_qsm '/inversion/sus_esharp_xy.nii']);
%     nii = make_nii(permute(sus,[1 3 2]),voxelSize);
%     save_nii(nii,[path_qsm '/inversion/sus_esharp_xz.nii']);
%     nii = make_nii(permute(sus,[2 3 1]),voxelSize);
%     save_nii(nii,[path_qsm '/inversion/sus_esharp_yz.nii']);
% end




%% save all variables for debugging purpose
if save_all
    clear nii;
    save('all.mat','-v7.3');
end

% save parameters used in the recon
save('parameters.mat','options','-v7.3')


%% clean up
% unix('rm *.nii*');
cd(init_dir);

