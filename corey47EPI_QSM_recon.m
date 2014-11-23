% 4.7T corey's EPI , 333 reosltuion

opt=se_epi_dw_optset();
opt.multivol_flg=1;
opt.upsmple_flag=0;
%opt.nii_z_invert=1;
opt.noask = 1; % automatically overwrite existing files (1), automatically create copies (2), or ask for each case (0)
opt.save_type=0;
opt.motion_c=0;
opt.only_basic = 0
opt.rcvr_comb=0
opt.homod = 0;


[par,img_out] = se_epi_dw_recon(pwd,opt);
matlabpool close;


[nPE, nRO, nSL, ne, nrcvrs] = size(img_out);

% voxelSize = [1.5 1.5 1.5];
% voxelSize = [2 2 2];
voxelSize = [par.lpe/nPE*10, par.lro/nRO*10, par.thk];

% !!!!!!!!!! flip into 4.7T scanner frame (coordinates)
img_all = flipdim(flipdim(img_out,1),2);



for i = 1:size(img_all,4) % all 320 time series

	img = squeeze(img_all(:,:,:,i,:));

	% matlabpool close;
	img_cmb = sense_se(img,voxelSize,[],3,3);

	img = img_cmb;
	mkdir('combine');
	nii = make_nii(abs(img),voxelSize);
	save_nii(nii,['combine/mag_cmb' num2str(i,'%03i') '.nii']);
	nii = make_nii(angle(img),voxelSize);
	save_nii(nii,['combine/ph_cmb' num2str(i,'%03i') '.nii']);


	bet_thr = 0.3;
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


	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% % unwrap combined phase with PRELUDE
	% disp('--> unwrap aliasing phase ...');
	% unix('prelude -a combine/mag_cmb${time_series}.nii -p combine/ph_cmb${time_series}.nii -u unph${time_series}.nii -m BET${time_series}_mask.nii -n 8');
	% unix('gunzip -f unph${time_series}.nii.gz');
	% nii = load_nii(['unph' num2str(i,'%03i') '.nii']);
	% unph = double(nii.img);
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



	smv_rad = 6;
	tik_reg = 0.001;
	% background field removal
	disp('--> RESHARP to remove background field ...');
	mkdir('RESHARP');
	[lph_resharp,mask_resharp] = resharp(unph,mask,voxelSize,smv_rad,tik_reg);
	lfs_resharp = -lph_resharp/(2.675e8*par.te*4.7)*1e6;
	nii = make_nii(lfs_resharp,voxelSize);
	save_nii(nii,['RESHARP/lfs_resharp' num2str(i,'%03i') '.nii']);

	lfs_poly = zeros(size(lfs_resharp));
	for sl = 1:size(lfs_resharp,3)
		lfs_poly(:,:,sl) = poly2d(lfs_resharp(:,:,sl),mask_resharp(:,:,sl));
	end

	nii = make_nii(lfs_poly,voxelSize);
	save_nii(nii,['RESHARP/lfs_resharp_poly' num2str(i,'%03i') '.nii']);


%% intrinsic euler angles 
% z-x-z convention, psi first, then theta, lastly phi
% psi and theta are left-handed, while gamma is right-handed!
	alpha = - par.psi/180*pi;
	beta = - par.theta/180*pi;
	gamma =  par.phi/180*pi;
	nor_vec = [sin(beta)*sin(gamma), sin(beta)*cos(gamma), cos(beta)]
	if ~ isequal(nor_vec,[0 0 1])
		disp('This is angled slicing');
		pwd
	end

	tvdi_n = 200;
	tv_reg = 3e-4;

	disp('--> TV susceptibility inversion ...');
	sus_resharp = tvdi(lfs_poly,mask_resharp,voxelSize,tv_reg,abs(img),nor_vec,tvdi_n);
	nii = make_nii(sus_resharp.*mask_resharp,voxelSize);
	save_nii(nii,['RESHARP/sus_resharp' num2str(i,'%03i') '.nii']);


end
