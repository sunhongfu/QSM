function qsm_swi47(path_in, path_out, options)
%QSM_SWI47 Quantitative susceptibility mapping from SWI sequence at 4.7T.
%   QSM_SWI47(PATH_IN, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH_IN     - directory of .fid from ge3d sequence      : ge3d__01.fid
%   PATH_OUT    - directory to save nifti and/or matrixes   : QSM_SWI_vxxx
%   OPTIONS     - parameter structure including fields below
%    .ref_coil  - reference coil to use for phase combine   : 3
%    .eig_rad   - radius (mm) of eig decomp kernel          : 3
%    .ph_unwrap - 'prelude' or 'laplacian' or 'bestpath'    : 'prelude'
%    .bkg_rm    - background field removal method(s)        : 'resharp'
%    .smv_rad   - radius (mm) of SMV convolution kernel     : 4
%    .tik_reg   - Tikhonov regularization for resharp       : 0.0005
%    .tv_reg    - Total variation regularization parameter  : 0.0005
%    .bet_thr   - threshold for BET brain mask              : 0.3
%    .tvdi_n    - iteration number of TVDI (nlcg)           : 200
%    .save_all   - save all the variables for debug          : 1


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
    options.ref_coil = 3;
end

if ~ isfield(options,'eig_rad')
    options.eig_rad = 4;
end

if ~ isfield(options,'bet_thr')
    options.bet_thr = 0.3;
end

if ~ isfield(options,'ph_unwrap')
    options.ph_unwrap = 'prelude';
    % % another option is
    % options.ph_unwrap = 'laplacian';
    % % prelude is preferred, unless there's sigularities
    % % in that case, have to use laplacian
    % another option is 'bestpath'
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
    options.inv_num = 200;
end

if ~ isfield(options,'save_all')
    options.save_all = 1;
end

if ~ isfield(options,'swi_ver')
    options.swi_ver = 'amir';
end

ref_coil  = options.ref_coil;
eig_rad   = options.eig_rad;
bet_thr   = options.bet_thr;
ph_unwrap = options.ph_unwrap;
bkg_rm    = options.bkg_rm;
smv_rad   = options.smv_rad;
tik_reg   = options.tik_reg;
t_svd     = options.t_svd;
tv_reg    = options.tv_reg;
inv_num   = options.inv_num;
save_all  = options.save_all;
swi_ver   = options.swi_ver;


%%% define directories
if strcmpi(ph_unwrap,'prelude')
    path_qsm = [path_out '/QSM_SWI47_v500'];
elseif strcmpi(ph_unwrap,'laplacian')
    path_qsm = [path_out '/QSM_SWI47_v500_lap'];
elseif strcmpi(ph_unwrap,'bestpath')
    path_qsm = [path_out '/QSM_SWI47_v500_best'];
end
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);


%%% generate raw img
disp('--> reconstruct fid to complex img ...');
[img,Pars] = swi47_recon(path_fid,swi_ver);


% scanner frame
img = permute(img, [2 1 3 4]);
img = flipdim(flipdim(img,2),3);
[nv,np,ns,~] = size(img); % phase, readout, slice, receivers
voxelSize = [Pars.lpe/nv, Pars.lro/np, Pars.lpe2/ns]*10;

% field directions
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


% combine receivers
if Pars.RCVRS_ > 1
    % combine RF coils
    disp('--> combine RF rcvrs ...');
    img_cmb = coils_cmb(img,voxelSize,ref_coil,eig_rad);
else  % single channel  
    img_cmb = img;
end

% save nifti
mkdir('combine');
nii = make_nii(abs(img_cmb),voxelSize);
save_nii(nii,'combine/mag_cmb.nii');
nii = make_nii(angle(img_cmb),voxelSize);
save_nii(nii,'combine/ph_cmb.nii');

clear img;


% generate brain mask
disp('--> extract brain volume and generate mask ...');
setenv('bet_thr',num2str(bet_thr));
unix('bet combine/mag_cmb.nii BET -f ${bet_thr} -m -R');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);


% phase unwrapping, prelude is preferred!
if strcmpi('prelude',ph_unwrap)
    % unwrap combined phase with PRELUDE
    disp('--> unwrap aliasing phase ...');
    bash_script = ['prelude -a combine/mag_cmb.nii -p combine/ph_cmb.nii ' ...
        '-u unph.nii -m BET_mask.nii -n 8'];
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
    Options.voxelSize = voxelSize;
    unph = lapunwrap(angle(img_cmb), Options);
    nii = make_nii(unph, voxelSize);
    save_nii(nii,'unph_lap.nii');

elseif strcmpi('bestpath',ph_unwrap)
    % unwrap the phase using best path
    fid = fopen('wrapped_phase.dat','w');
    fwrite(fid,angle(img_cmb),'float');
    fclose(fid);
    % mask_unwrp = uint8(hemo_mask.*255);
    mask_unwrp = uint8(abs(mask)*255);
    fid = fopen('mask_unwrp.dat','w');
    fwrite(fid,mask_unwrp,'uchar');
    fclose(fid);

    unix('cp /home/hongfu/Documents/MATLAB/3DSRNCP 3DSRNCP');
    setenv('nv',num2str(nv));
    setenv('np',num2str(np));
    setenv('ns',num2str(ns));
    bash_script = ['./3DSRNCP wrapped_phase.dat mask_unwrp.dat unwrapped_phase.dat ' ...
        '$nv $np $ns reliability.dat'];
    unix(bash_script) ;

    fid = fopen('unwrapped_phase.dat','r');
    unph = fread(fid,'float');
    unph = reshape(unph - unph(1) ,[nv,np,ns]);
    fclose(fid);
    nii = make_nii(unph,voxelSize);
    save_nii(nii,'unph.nii');

    % fid = fopen('reliability.dat','r');
    % reliability = fread(fid,'float');
    % fclose(fid);
    % reliability = reshape(reliability,[nv,np,ns]);
    % reliability = 1./reliability.*mask;
    % reliability_smooth = smooth3(reliability,'box',[7,7,3]);
    % % reliability(reliability <= 0.05) = 0;
    % % reliability(reliability > 0.05) = 1;
    % nii = make_nii(reliability_smooth,voxelSize);
    % save_nii(nii,'reliability_smooth.nii');
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


% PDF
if sum(strcmpi('pdf',bkg_rm))
    disp('--> PDF to remove background field ...');
    [lfs_pdf,mask_pdf] = pdf(tfs,mask,voxelSize,smv_rad, ...
        abs(img_cmb),z_prjs);
    % 3D 2nd order polyfit to remove any residual background
    lfs_pdf= poly3d(lfs_pdf,mask_pdf);

    % save nifti
    mkdir('PDF');
    nii = make_nii(lfs_pdf,voxelSize);
    save_nii(nii,'PDF/lfs_pdf.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on PDF...');
    sus_pdf = tvdi(lfs_pdf, mask_pdf, voxelSize, tv_reg, ...
        abs(img_cmb), z_prjs, inv_num); 

    % save nifti
    nii = make_nii(sus_pdf.*mask_pdf,voxelSize);
    save_nii(nii,'PDF/sus_pdf.nii');
end


% SHARP (t_svd: truncation threthold for t_svd)
if sum(strcmpi('sharp',bkg_rm))
    disp('--> SHARP to remove background field ...');
    [lfs_sharp, mask_sharp] = sharp(tfs,mask,voxelSize,smv_rad,t_svd);
    % 3D 2nd order polyfit to remove any residual background
    lfs_sharp= poly3d(lfs_sharp,mask_sharp);

    % save nifti
    mkdir('SHARP');
    nii = make_nii(lfs_sharp,voxelSize);
    save_nii(nii,'SHARP/lfs_sharp.nii');
    
    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on SHARP...');
    sus_sharp = tvdi(lfs_sharp, mask_sharp, voxelSize, tv_reg, ...
        abs(img_cmb), z_prjs, inv_num); 
   
    % save nifti
    nii = make_nii(sus_sharp.*mask_sharp,voxelSize);
    save_nii(nii,'SHARP/sus_sharp.nii');
end


% RE-SHARP (tik_reg: Tikhonov regularization parameter)
if sum(strcmpi('resharp',bkg_rm))
    disp('--> RESHARP to remove background field ...');
    [lfs_resharp, mask_resharp] = resharp(tfs,mask,voxelSize,smv_rad,tik_reg);
    % 3D 2nd order polyfit to remove any residual background
    lfs_resharp= poly3d(lfs_resharp,mask_resharp);

    % save nifti
    mkdir('RESHARP');
    nii = make_nii(lfs_resharp,voxelSize);
    save_nii(nii,'RESHARP/lfs_resharp.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on RESHARP...');
    sus_resharp = tvdi(lfs_resharp, mask_resharp, voxelSize, tv_reg, ...
        abs(img_cmb), z_prjs, inv_num); 
   
    % save nifti
    nii = make_nii(sus_resharp.*mask_resharp,voxelSize);
    save_nii(nii,'RESHARP/sus_resharp.nii');

end


% LBV
if sum(strcmpi('lbv',bkg_rm))
    disp('--> LBV to remove background field ...');
    lfs_lbv = LBV(tfs,mask,size(tfs),voxelSize,0.01,2); % strip 2 layers
    mask_lbv = ones(size(mask));
    mask_lbv(lfs_lbv==0) = 0;
    % 3D 2nd order polyfit to remove any residual background
    lfs_lbv= poly3d(lfs_lbv,mask_lbv);

    % save nifti
    mkdir('LBV');
    nii = make_nii(lfs_lbv,voxelSize);
    save_nii(nii,'LBV/lfs_lbv.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on lbv...');
    sus_lbv = tvdi(lfs_lbv,mask_lbv,voxelSize,tv_reg, ...
        abs(img_cmb),z_prjs,inv_num);   

    % save nifti
    nii = make_nii(sus_lbv.*mask_lbv,voxelSize);
    save_nii(nii,'LBV/sus_lbv.nii');

end


% save all variables for debugging purpose
if save_all
    clear nii;
    save('all.mat','-v7.3');
end

% save parameters used in the recon
save('parameters.mat','options','-v7.3')


% clean up
% unix('rm *.nii*');
cd(init_dir);
