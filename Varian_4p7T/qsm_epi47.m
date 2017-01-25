function qsm_epi47(path_in, path_out, options)
%QSM_EPI47 Quantitative susceptibility mapping from EPI (Corey's) at 4.7T.
%   QSM_EPI47(PATH_IN, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH_IN      - directory of .fid from gemsme3d sequence  : se_epi_dw***.fid
%   PATH_OUT     - directory to save nifti and/or matrixes   : QSM_EPI
%   OPTIONS      - parameter structure including fields below
%    .ref_coil   - reference coil to use for phase combine   : 3
%    .eig_rad    - radius (mm) of eig decomp kernel          : 5
%    .bet_thr    - threshold for BET brain mask              : 0.35
%    .bet_smooth - smoothness of BET brain mask at edges     : 2
%    .ph_unwrap  - 'prelude' or 'laplacian' or 'bestpath'    : 'prelude'
%    .bkg_rm     - background field removal method(s)        : 'resharp'
%    .smv_rad    - radius (mm) of SMV convolution kernel     : 5
%    .tik_reg    - Tikhonov regularization for RESHARP       : 5e-4
%    .t_svd      - truncation of SVD for SHARP               : 0.05
%    .lbv_layer  - number of layers to be stripped off LBV   : 1
%    .tv_reg     - Total variation regularization parameter  : 5e-4
%    .inv_num    - iteration number of TVDI (nlcg)           : 500
%    .save_all   - save the entire workspace/variables       : 0

% default settings
if ~ exist('path_in','var') || isempty(path_in)
    path_in = pwd;
end

if exist([path_in '/fid'],'file')
    path_fid = path_in;
    path_fid = cd(cd(path_fid));
elseif exist(path_in,'file')
	[pathstr,name,ext] = fileparts(path_in);
	path_fid = cd(cd(pathstr));
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
    options.eig_rad = 5;
end

if ~ isfield(options,'bet_thr')
    options.bet_thr = 0.35;
end

if ~ isfield(options,'bet_smooth')
    options.bet_smooth = 2;
end

if ~ isfield(options,'ph_unwrap')
    % options.ph_unwrap = 'prelude';
    % % another option is
    options.ph_unwrap = 'laplacian';
    % % prelude is preferred, unless there's sigularities
    % % in that case, have to use laplacian
    % another option is 'bestpath'
end

if ~ isfield(options,'bkg_rm')
    options.bkg_rm = 'resharp';
    % options.bkg_rm = {'resharp','lbv'};
    % options.bkg_rm = {'pdf','sharp','resharp','lbv'};
end

if ~ isfield(options,'t_svd')
    options.t_svd = 0.05;
end

if ~ isfield(options,'smv_rad')
    options.smv_rad = 5;
end

if ~ isfield(options,'tik_reg')
    options.tik_reg = 5e-4;
end

if ~ isfield(options,'lbv_layer')
    options.lbv_layer = 1;
end

if ~ isfield(options,'tv_reg')
    options.tv_reg = 5e-4;
end

if ~ isfield(options,'inv_num')
    options.inv_num = 500;
end

if ~ isfield(options,'save_all')
    options.save_all = 0;
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
save_all   = options.save_all;


% 4.7T corey's EPI parameters
opt=se_epi_dw_optset();
opt.multivol_flg=1;
opt.upsmple_flag=0;
opt.noask = 1; 
opt.save_type=0;
opt.motion_c=0;
opt.only_basic = 0;
opt.rcvr_comb=0;
opt.homod = 0;


%% define directories
if strcmpi(ph_unwrap,'prelude')
    path_qsm = [path_out '/QSM_EPI47_pre'];
elseif strcmpi(ph_unwrap,'laplacian')
    path_qsm = [path_out '/QSM_EPI47_lap'];
elseif strcmpi(ph_unwrap,'bestpath')
    path_qsm = [path_out '/QSM_EPI47_best'];
end
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);


% reconstruct complex image from fid file
[par,img_out] = se_epi_dw_recon(path_fid,opt);
matlabpool close


% flip to match 4.7T scanner frame/gradients (coordinates)
[nv, np, ns, ~, ~] = size(img_out);
voxelSize = [par.lpe/nv*10, par.lro/np*10, par.thk];
img_all = flipdim(flipdim(img_out,1),2); % the same as rot180 (rot90(x,2))


% a quick peak at the raw phase
disp('check out the uncombined raw phase!');
nii = make_nii(squeeze(angle(img_all(:,:,:,1,:))),voxelSize);
save_nii(nii,'rawphase.nii');

[nv np ns nr nrcvrs] = size(img_all);
% nr is the number of runs

% recon combined magnitude
mkdir('combine');
disp('--> combine multiple channels ...');
for i = 1:nr % all time series
	img = squeeze(img_all(:,:,:,i,:));
    
	if par.nrcvrs > 1
		img_cmb = adaptive_cmb(img,voxelSize,ref_coil,eig_rad,0);
	else
		img_cmb = img;
	end

	img_cmb_all(:,:,:,i) = img_cmb;

	nii = make_nii(abs(img_cmb),voxelSize);
	save_nii(nii,['combine/mag_cmb' num2str(i,'%03i') '.nii']);
	nii = make_nii(angle(img_cmb),voxelSize);
	save_nii(nii,['combine/ph_cmb' num2str(i,'%03i') '.nii']);
end
nii = make_nii(abs(img_cmb_all),voxelSize);
save_nii(nii,'all_mag_cmb.nii');

% generate BET on 1st volume
[status,cmdout] = unix('rm BET*');
disp('--> extract brain volume and generate mask ...');
setenv('bet_thr',num2str(bet_thr));
setenv('bet_smooth',num2str(bet_smooth));
bash_script = ['bet2 combine/mag_cmb001.nii BET ' ...
	'-f ${bet_thr} -m'];
unix(bash_script);
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii(['BET_mask.nii']);
mask = double(nii.img);

mask_rep = repmat(mask,[1 1 1 nr]);
nii = make_nii(mask_rep,voxelSize);
save_nii(nii,'mask_rep.nii');

% spm to align all volumes
P = spm_select('ExtList', pwd, '^all_mag_cmb.nii',Inf);
flags.mask=0;
spm_realign(P);
load all_mag_cmb.mat
m=[voxelSize(1) 0 0 0; 0 voxelSize(2) 0 0; 0 0 voxelSize(3) 0; 0 0 0 1];
for i = 2:size(mat,3)
	% mat(:,:,i) = inv(inv(m)*mat(:,:,i))*m;
		mat(:,:,i) = m*inv(mat(:,:,i))*m;
end
save('mask_rep.mat','mat');
P = spm_select('ExtList', pwd, '^mask_rep.nii',Inf);
flags.mask=0;
spm_reslice(P,flags);

nii = load_nii('rmask_rep.nii');
rmask = nii.img;
rmask(isnan(rmask)) = 0;
rmask(isinf(rmask)) = 0;

% process QSM on individual run volume
for i = 1:nr % all time series
	img_cmb = img_cmb_all(:,:,:,i);
	nii = make_nii(rmask(:,:,:,i),voxelSize);
	save_nii(nii,['BET' num2str(i,'%03i') '_mask.nii']);

	% unwrap the phase
	if strcmpi('prelude',ph_unwrap)
	    % unwrap combined phase with PRELUDE
	    disp('--> unwrap aliasing phase ...');
	    setenv('time_series',num2str(i,'%03i'));
	    bash_script = ['prelude -a combine/mag_cmb${time_series}.nii ' ...
	        '-p combine/ph_cmb${time_series}.nii -u unph${time_series}.nii ' ...
	        '-m BET${time_series}_mask.nii -n 8'];
	    unix(bash_script);
	    unix('gunzip -f unph${time_series}.nii.gz');
	    nii = load_nii(['unph' num2str(i,'%03i') '.nii']);
	    unph = double(nii.img);

	elseif strcmpi('laplacian',ph_unwrap)
		% Ryan Topfer's Laplacian unwrapping
	    disp('--> unwrap aliasing phase using laplacian...');
	    Options.voxelSize = voxelSize;
		unph = lapunwrap(angle(img_cmb), Options);
		nii = make_nii(unph, voxelSize);
		save_nii(nii,['unph_lap' num2str(i,'%03i') '.nii']);

	elseif strcmpi('bestpath',ph_unwrap)
	    % unwrap the phase using best path
	    disp('--> unwrap aliasing phase using bestpath...');
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
	    unph = reshape(unph - unph(1) ,[nv, np, ns]);
	    fclose(fid);
	    nii = make_nii(unph,voxelSize);
	    save_nii(nii,'unph.nii');

	end


	% normalize to ppm unit
	tfs = -unph/(2.675e8*par.te*4.7)*1e6; % unit ppm
	nii = make_nii(tfs,voxelSize);
	save_nii(nii,'tfs.nii');


	% intrinsic euler angles 
	% z-x-z convention, psi first, then theta, lastly phi
	% psi and theta are left-handed, while gamma is right-handed!
	% alpha = - par.psi/180*pi;
	beta = - par.theta/180*pi;
	gamma =  par.phi/180*pi;
	z_prjs = [sin(beta)*sin(gamma), sin(beta)*cos(gamma), cos(beta)];
	if ~ isequal(z_prjs,[0 0 1])
		disp('This is angled slicing');
		disp(z_prjs);
	end


	% PDF
	if sum(strcmpi('pdf',bkg_rm))
	    disp('--> PDF to remove background field ...');
	    [lfs_pdf,mask_pdf] = pdf(tfs,mask,voxelSize,smv_rad, ...
	        abs(img_cmb),z_prjs);
	    % 2D 2nd order polyfit to remove any residual background
	    lfs_pdf = lfs_pdf - poly2d(lfs_pdf,mask_pdf);

	    % save nifti
	    mkdir('PDF');
	    nii = make_nii(lfs_pdf,voxelSize);
    	save_nii(nii,['PDF/lfs_pdf' num2str(i,'%03i') '.nii']);

	    % inversion of susceptibility 
	    disp('--> TV susceptibility inversion on PDF...');
	    sus_pdf = tvdi(lfs_pdf, mask_pdf, voxelSize, tv_reg, ...
	        abs(img_cmb), z_prjs, inv_num); 

	   	% demean susceptibility
	   	sus_pdf = sus_pdf - mean(sus_pdf(logical(mask_pdf(:))));

	    % save nifti
	    nii = make_nii(sus_pdf.*mask_pdf,voxelSize);
	    save_nii(nii,['PDF/sus_pdf' num2str(i,'%03i') '.nii']);

	    if save_all
			lfs_pdf_all(:,:,:,i) = lfs_pdf;
			mask_pdf_all(:,:,:,i) = mask_pdf;
			sus_pdf_all(:,:,:,i) = sus_pdf;
		end
	end


	% SHARP (t_svd: truncation threthold for t_svd)
	if sum(strcmpi('sharp',bkg_rm))
	    disp('--> SHARP to remove background field ...');
	    [lfs_sharp, mask_sharp] = sharp(tfs,mask,voxelSize,smv_rad,t_svd);
	    % 2D 2nd order polyfit to remove any residual background
	    lfs_sharp = lfs_sharp - poly2d(lfs_sharp,mask_sharp);

	    % save nifti
	    mkdir('SHARP');
	    nii = make_nii(lfs_sharp,voxelSize);
	    save_nii(nii,['SHARP/lfs_sharp' num2str(i,'%03i') '.nii']);
	    
	    % inversion of susceptibility 
	    disp('--> TV susceptibility inversion on SHARP...');
	    sus_sharp = tvdi(lfs_sharp, mask_sharp, voxelSize, tv_reg, ...
	        abs(img_cmb), z_prjs, inv_num); 

	   	% demean susceptibility
	   	sus_sharp = sus_sharp - mean(sus_sharp(logical(mask_sharp(:))));
	   
	    % save nifti
	    nii = make_nii(sus_sharp.*mask_sharp,voxelSize);
	    save_nii(nii,['SHARP/sus_sharp' num2str(i,'%03i') '.nii']);

	    if save_all
			lfs_sharp_all(:,:,:,i) = lfs_sharp;
			mask_sharp_all(:,:,:,i) = mask_sharp;
			sus_sharp_all(:,:,:,i) = sus_sharp;
		end
	end


	% RE-SHARP (tik_reg: Tikhonov regularization parameter)
	if sum(strcmpi('resharp',bkg_rm))
	    disp('--> RESHARP to remove background field ...');
	    [lfs_resharp, mask_resharp] = resharp(tfs,mask,voxelSize,smv_rad,tik_reg);
	    % 2D 2nd order polyfit to remove any residual background
	    lfs_resharp = lfs_resharp - poly2d(lfs_resharp,mask_resharp);

	    % save nifti
	    mkdir('RESHARP');
	    nii = make_nii(lfs_resharp,voxelSize);
	    save_nii(nii,['RESHARP/lfs_resharp' num2str(i,'%03i') '.nii']);

	    % inversion of susceptibility 
	    disp('--> TV susceptibility inversion on RESHARP...');
	    sus_resharp = tvdi(lfs_resharp, mask_resharp, voxelSize, tv_reg, ...
	        abs(img_cmb), z_prjs, inv_num); 
	   	
	   	% demean susceptibility
	   	sus_resharp = sus_resharp - mean(sus_resharp(logical(mask_resharp(:))));

	    % save nifti
	    nii = make_nii(sus_resharp.*mask_resharp,voxelSize);
	    save_nii(nii,['RESHARP/sus_resharp' num2str(i,'%03i') '.nii']);

	    if save_all
			lfs_resharp_all(:,:,:,i) = lfs_resharp;
			mask_resharp_all(:,:,:,i) = mask_resharp;
			sus_resharp_all(:,:,:,i) = sus_resharp;
		end
	end
	
	
	if sum(strcmpi('esharp',bkg_rm))
	    disp('--> E-SHARP to remove background field ...');
	    Parameters.voxelSize             = voxelSize; % in mm
	    Parameters.resharpRegularization = tik_reg ;
	    Parameters.resharpKernelRadius   = smv_rad ; % in mm
	    Parameters.radius                = [ 10 10 5 ] ;


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
	   
	    lfs_esharp      = localField;
	    mask_esharp     = mask;


	    % 2D 2nd order polyfit to remove any residual background
	    lfs_esharp = lfs_esharp - poly2d(localField,mask_esharp);

	    % save nifti
	    mkdir('ESHARP');
	    nii = make_nii(lfs_esharp,voxelSize);
	    save_nii(nii,['ESHARP/lfs_esharp' num2str(i,'%03i') '.nii']);

	    % inversion of susceptibility 
	    disp('--> TV susceptibility inversion on RESHARP...');
	    sus_esharp = tvdi(lfs_esharp, mask_esharp, voxelSize, tv_reg, ...
	        abs(img_cmb), z_prjs, inv_num); 
	   	
	   	% demean susceptibility
	   	sus_esharp = sus_esharp - mean(sus_esharp(logical(mask_esharp(:))));

	    % save nifti
	    nii = make_nii(sus_esharp.*mask_esharp,voxelSize);
	    save_nii(nii,['ESHARP/sus_esharp' num2str(i,'%03i') '.nii']);

	    if save_all
			lfs_esharp_all(:,:,:,i) = lfs_esharp;
			mask_esharp_all(:,:,:,i) = mask_esharp;
			sus_esharp_all(:,:,:,i) = sus_esharp;
		end
	end


	% LBV
	if sum(strcmpi('lbv',bkg_rm))
	    disp('--> LBV to remove background field ...');
	    lfs_lbv = LBV(tfs,mask,size(tfs),voxelSize,0.01,lbv_layer); % strip 2 layers
	    mask_lbv = ones(size(mask));
	    mask_lbv(lfs_lbv==0) = 0;
	    % 2D 2nd order polyfit to remove any residual background
	    lfs_lbv = lfs_lbv - poly2d(lfs_lbv,mask_lbv);

	    % save nifti
	    mkdir('LBV');
	    nii = make_nii(lfs_lbv,voxelSize);
	    save_nii(nii,['LBV/lfs_lbv' num2str(i,'%03i') '.nii']);

	    % inversion of susceptibility 
	    disp('--> TV susceptibility inversion on lbv...');
	    sus_lbv = tvdi(lfs_lbv,mask_lbv,voxelSize,tv_reg, ...
	        abs(img_cmb),z_prjs,inv_num);   

	   	% demean susceptibility
	   	sus_lbv = sus_lbv - mean(sus_lbv(logical(mask_lbv(:))));

	    % save nifti
	    nii = make_nii(sus_lbv.*mask_lbv,voxelSize);
	    save_nii(nii,['LBV/sus_lbv' num2str(i,'%03i') '.nii']);

	    if save_all
			lfs_lbv_all(:,:,:,i) = lfs_lbv;
			mask_lbv_all(:,:,:,i) = mask_lbv;
			sus_lbv_all(:,:,:,i) = sus_lbv;
		end
	end


	% to save all the variables
	if save_all
		mask_all(:,:,:,i) = mask;
		unph_all(:,:,:,i) = unph;
	end

end


% save all variables for debugging purpose
if save_all
    clear nii;
    save('all.mat','-v7.3');
end

% save parameters used in the recon
save('parameters.mat','options','-v7.3')

% save the git log for future tracking
unix('git log --branches --decorate --color --abbrev-commit --graph --no-merges --tags > git_log');

% clean up
unix('rm *.nii*');
cd(init_dir);
