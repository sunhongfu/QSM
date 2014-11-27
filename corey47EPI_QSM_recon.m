function qsm_epi47(path_in, path_out, options)
%QSM_EPI47 Quantitative susceptibility mapping from EPI sequence (Corey's) at 4.7T.
%   QSM_EPI47(PATH_IN, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH_IN    - directory of .fid from gemsme3d sequence  : se_epi_dw***.fid
%   PATH_OUT   - directory to save nifti and/or matrixes   : QSM_EPI_vxxx
%   OPTIONS    - parameter structure including fields below
%    .bet_thr  - threshold for BET brain mask              : 0.3
%    .bkg_rm   - background field removal method(s)        : 'resharp'
%    .smv_rad  - radius (mm) of SMV convolution kernel     : 4
%    .tik_reg  - Tikhonov regularization for RESHARP       : 0.0005
%    .t_svd    - truncation of SVD for SHARP               : 0.05
%    .tv_reg   - Total variation regularization parameter  : 0.0005
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
    options.ref_coil = 3;
end

if ~ isfield(options,'eig_rad')
    options.eig_rad = 4;
end

if ~ isfield(options,'bet_thr')
    options.bet_thr = 0.3;
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


ref_coil = options.ref_coil;
eig_rad  = options.eig_rad;
bet_thr  = options.bet_thr;
bkg_rm   = options.bkg_rm;
smv_rad  = options.smv_rad;
tik_reg  = options.tik_reg;
% t_svd    = options.t_svd;
tv_reg   = options.tv_reg;
inv_num  = options.inv_num;
save_all = options.save_all;


% 4.7T corey's EPI parameters
opt=se_epi_dw_optset();
opt.multivol_flg=1;
opt.upsmple_flag=0;
opt.noask = 1; 
opt.save_type=0;
opt.motion_c=0;
opt.only_basic = 0
opt.rcvr_comb=0
opt.homod = 0;


% define directories
path_qsm = [path_out '/QSM_EPI_v100'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);



% reconstruct complex image from fid file
[par,img_out] = se_epi_dw_recon(pwd,opt);
matlabpool close;


% !!!!!!!!!! flip to match 4.7T scanner frame (coordinates)
[nPE, nRO, nSL, nRN, nrcvrs] = size(img_out);
voxelSize = [par.lpe/nPE*10, par.lro/nRO*10, par.thk];
img_all = flipdim(flipdim(img_out,1),2);


% process QSM on individual run volume
for i = 1:size(img_all,4) % all time series
	img = squeeze(img_all(:,:,:,i,:));

	if par.nrcvrs > 1
		img = coils_cmb(img,voxelSize,ref_coil,eig_rad);
	end

	mkdir('combine');
	nii = make_nii(abs(img),voxelSize);
	save_nii(nii,['combine/mag_cmb' num2str(i,'%03i') '.nii']);
	nii = make_nii(angle(img),voxelSize);
	save_nii(nii,['combine/ph_cmb' num2str(i,'%03i') '.nii']);


	disp('--> extract brain volume and generate mask ...');
	setenv('bet_thr',num2str(bet_thr));
	setenv('time_series',num2str(i,'%03i'));
	unix('bet combine/mag_cmb${time_series}.nii BET${time_series} -f ${bet_thr} -m -Z');
	unix('gunzip -f BET${time_series}.nii.gz');
	unix('gunzip -f BET${time_series}_mask.nii.gz');
	nii = load_nii(['BET' num2str(i,'%03i') '_mask.nii']);
	mask = double(nii.img);


	Options.voxelSize = voxelSize;
	unph = lapunwrap(angle(img), Options);
	nii = make_nii(unph, voxelSize);
	save_nii(nii,['unph_lap' num2str(i,'%03i') '.nii']);


	% background field removal
	disp('--> RESHARP to remove background field ...');
	mkdir('RESHARP');
	[lph_resharp,mask_resharp] = resharp(unph,mask,voxelSize,smv_rad,tik_reg);
	lfs_resharp = -lph_resharp/(2.675e8*par.te*4.7)*1e6;
	% nii = make_nii(lfs_resharp,voxelSize);
	% save_nii(nii,['RESHARP/lfs_resharp' num2str(i,'%03i') '.nii']);

	lfs_poly = zeros(size(lfs_resharp));
	for sl = 1:size(lfs_resharp,3)
		lfs_poly(:,:,sl) = poly2d(lfs_resharp(:,:,sl),mask_resharp(:,:,sl));
	end

	nii = make_nii(lfs_poly,voxelSize);
	save_nii(nii,['RESHARP/lfs_resharp_poly' num2str(i,'%03i') '.nii']);


	% intrinsic euler angles 
	% z-x-z convention, psi first, then theta, lastly phi
	% psi and theta are left-handed, while gamma is right-handed!
	% alpha = - par.psi/180*pi;
	beta = - par.theta/180*pi;
	gamma =  par.phi/180*pi;
	z_prjs = [sin(beta)*sin(gamma), sin(beta)*cos(gamma), cos(beta)]
	if ~ isequal(z_prjs,[0 0 1])
		disp('This is angled slicing');
		disp(z_prjs);
	end

	disp('--> TV susceptibility inversion ...');
	sus_resharp = tvdi(lfs_poly,mask_resharp,voxelSize,tv_reg,abs(img),z_prjs,inv_num);
	nii = make_nii(sus_resharp.*mask_resharp,voxelSize);
	save_nii(nii,['RESHARP/sus_resharp' num2str(i,'%03i') '.nii']);


end
