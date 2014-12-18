function qsm_epi47(path_in, path_out, options)
%QSM_EPI47 Quantitative susceptibility mapping from EPI (Corey's) at 4.7T.
%   QSM_EPI47(PATH_IN, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH_IN    - directory of .fid from gemsme3d sequence  : se_epi_dw***.fid
%   PATH_OUT   - directory to save nifti and/or matrixes   : QSM_EPI_vxxx
%   OPTIONS    - parameter structure including fields below
%    .bet_thr  - threshold for BET brain mask              : 0.3
%    .smv_rad  - radius (mm) of SMV convolution kernel     : 4
%    .tik_reg  - Tikhonov regularization for RESHARP       : 5e-4
%    .tv_reg   - Total variation regularization parameter  : 5e-4
%    .inv_num  - iteration number of TVDI (nlcg)           : 200
%    .save_all - save all the variables for debug (~ 0)    : 1

% default settings
if ~ exist('path_in','var') || isempty(path_in)
    path_in = pwd;
end

if exist([path_in '/fid'],'file')
    path_fid = path_in;
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
    options.ref_coil = 1;
end

if ~ isfield(options,'eig_rad')
    options.eig_rad = 4;
end

if ~ isfield(options,'bet_thr')
    options.bet_thr = 0.3;
end

if ~ isfield(options,'ph_unwrap')
    options.ph_unwrap = 'laplacian';
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

% if ~ isfield(options,'t_svd')
%     options.t_svd = 0.05;
% end

if ~ isfield(options,'tv_reg')
    options.tv_reg = 5e-4;
end

if ~ isfield(options,'inv_num')
    options.inv_num = 200;
end

if ~ isfield(options,'save_all')
    options.save_all = 1;
end


ref_coil  = options.ref_coil;
eig_rad   = options.eig_rad;
bet_thr   = options.bet_thr;
ph_unwrap = options.ph_unwrap;
% bkg_rm    = options.bkg_rm;
smv_rad   = options.smv_rad;
tik_reg   = options.tik_reg;
% t_svd     = options.t_svd;
tv_reg    = options.tv_reg;
inv_num   = options.inv_num;
save_all  = options.save_all;


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


%%% define directories
if strcmpi(ph_unwrap,'prelude')
    path_qsm = [path_out '/QSM_EPI15_v1'];
elseif strcmpi(ph_unwrap,'laplacian')
    path_qsm = [path_out '/QSM_EPI15_v1_lap'];
elseif strcmpi(ph_unwrap,'bestpath')
    path_qsm = [path_out '/QSM_EPI15_v1_best'];
end
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);


% reconstruct complex image from fid file
[par,img_out] = se_epi_dw_recon(path_fid,opt);
matlabpool close;


% flip to match 4.7T scanner frame/gradients (coordinates)
[nv, np, ns, nrcvrs, ~] = size(img_out);
voxelSize = [par.lpe/nv*10, par.lro/np*10, par.thk];
img_all = flipdim(flipdim(img_out,1),2); % the same as rot180 (rot90(x,2))

% save all the variables
if save_all
	img_cmb_all = zeros([nv, np, ns, nrcvrs]);
	mask_all = zeros([nv, np, ns, nrcvrs]);
	unph_all = zeros([nv, np, ns, nrcvrs]);
	lfs_resharp_all = zeros([nv, np, ns, nrcvrs]);
	mask_resharp_all = zeros([nv, np, ns, nrcvrs]);
	lfs_poly_all = zeros([nv, np, ns, nrcvrs]);
	sus_resharp_all = zeros([nv, np, ns, nrcvrs]);
end
    
% process QSM on individual run volume
for i = 1:size(img_all,4) % all time series
	img = squeeze(img_all(:,:,:,i,:));
    
    disp('--> combine multiple channels ...');
	if par.nrcvrs > 1
		img_cmb = coils_cmb(img,voxelSize,ref_coil,eig_rad);
	end

	mkdir('combine');
	nii = make_nii(abs(img_cmb),voxelSize);
	save_nii(nii,['combine/mag_cmb' num2str(i,'%03i') '.nii']);
	nii = make_nii(angle(img_cmb),voxelSize);
	save_nii(nii,['combine/ph_cmb' num2str(i,'%03i') '.nii']);

	disp('--> extract brain volume and generate mask ...');
	setenv('bet_thr',num2str(bet_thr));
	setenv('time_series',num2str(i,'%03i'));
	[status,cmdout] = unix('rm BET*');
	bash_script = ['bet combine/mag_cmb${time_series}.nii BET${time_series} ' ...
		'-f ${bet_thr} -m -Z'];
	unix(bash_script);
	unix('gunzip -f BET${time_series}.nii.gz');
	unix('gunzip -f BET${time_series}_mask.nii.gz');
	nii = load_nii(['BET' num2str(i,'%03i') '_mask.nii']);
	mask = double(nii.img);

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

	% background field removal
	disp('--> RESHARP to remove background field ...');
	mkdir('RESHARP');
	[lph_resharp,mask_resharp] = resharp(unph,mask,voxelSize,smv_rad,tik_reg);
	lfs_resharp = -lph_resharp/(2.675e8*par.te*4.7)*1e6;
	% nii = make_nii(lfs_resharp,voxelSize);
	% save_nii(nii,['RESHARP/lfs_resharp' num2str(i,'%03i') '.nii']);

	lfs_poly = poly2d(lfs_resharp,mask_resharp);
	nii = make_nii(lfs_poly,voxelSize);
	save_nii(nii,['RESHARP/lfs_resharp_poly' num2str(i,'%03i') '.nii']);

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

	disp('--> TV susceptibility inversion ...');
	sus_resharp = tvdi(lfs_poly,mask_resharp,voxelSize,tv_reg, ...
		abs(img_cmb),z_prjs,inv_num);
	nii = make_nii(sus_resharp.*mask_resharp,voxelSize);
	save_nii(nii,['RESHARP/sus_resharp' num2str(i,'%03i') '.nii']);

	% to save all the variables
	if save_all
		img_cmb_all(:,:,:,i) = img_cmb;
		mask_all(:,:,:,i) = mask;
		unph_all(:,:,:,i) = unph;
		lfs_resharp_all(:,:,:,i) = lfs_resharp;
		mask_resharp_all(:,:,:,i) = mask_resharp;
		lfs_poly_all(:,:,:,i) = lfs_poly;
		sus_resharp_all(:,:,:,i) = sus_resharp;
	end

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
