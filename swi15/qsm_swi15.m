function qsm_swi15(path_in)
% QSM reconstruction for 1.5T SWI
init_dir = pwd;

if ~ exist('path_in','var') || isempty(path_in)
	path_in = pwd;
end

listing = dir([path_in '/*.out']);
pathstr = cd(cd(path_in));

for i = 1:size(listing,1)

	filename = listing(i).name;

	rawfile = {[pathstr,filesep],filename};

	path_out = [filename,'_QSM'];
	
	mkdir(pathstr,path_out);

	cd([pathstr,filesep,path_out]);


	% generate raw img
	disp('--> (1/7) reconstruct to complex img ...');
	[img,params] = swi15_recon(rawfile);


	% size and resolution
	[Npe,Nro,Ns,~] = size(img);
	FOV = params.protocol_header.sSliceArray.asSlice{1};
	voxelSize = [FOV.dPhaseFOV/Npe, FOV.dReadoutFOV/Nro,  FOV.dThickness/Ns];



	% combine RF coils
	disp('--> (2/7) combine RF rcvrs ...');
	cref = 8; % reference coil
	radi = 5; % kernel size

	img_cmb = sense_se(img,voxelSize,cref,radi);
	nii = make_nii(abs(img_cmb),voxelSize);
	save_nii(nii,'mag_cmb.nii');
	nii = make_nii(angle(img_cmb),voxelSize);
	save_nii(nii,'ph_cmb.nii');



	% % combine coils (2D SVD)
	% img_cmb = zeros([Npe,Nro,Ns]);
	% matlabpool open
	% parfor i = 1:Ns
	%     img_cmb(:,:,i) = coilCombinePar(img(:,:,i,:));
	% end
	% matlabpool close
	% nii = make_nii(abs(img_cmb),voxelSize);
	% save_nii(nii,'mag_cmb.nii');
	% nii = make_nii(angle(img_cmb),voxelSize);
	% save_nii(nii,'ph_cmb.nii');



	% generate brain mask
	disp('--> (3/7) extract brain volume and generate mask ...');

	! bet mag_cmb.nii BET -f 0.4 -m -R;
	! gunzip -f BET.nii.gz;
	! gunzip -f BET_mask.nii.gz;
	nii = load_nii('BET_mask.nii');
	mask = double(nii.img);



	% unwrap combined phase with PRELUDE in 3D
	disp('--> (4/7) unwrap aliasing phase ...');

	! prelude -a mag_cmb.nii -p ph_cmb.nii -u unph.nii -m BET_mask.nii -n 8
	! gunzip -f unph.nii.gz
	nii = load_nii('unph.nii');
	unph = double(nii.img);


	% Options.voxelSize = voxelSize;
	% unph = lapunwrap(angle(img_cmb), Options);
	% nii = make_nii(unph,voxelSize);
	% save_nii(nii,'unph_lap.nii');


	% background field removal
	disp('--> (6/7) RESHARP to remove background field ...');
	ker_rad = 5; % convolution kernel radius size (mm)
	tik_reg = 1e-3; % tikhonov regularization

	% % (1) PDF
	% theta = -acos(params.protocol_header.sSliceArray.asSlice{1}.sNormal.dTra);
	% [lfs,mask_ero] = pdf(tfs,mask,voxelSize,ker_rad,abs(img_cmb),theta);
	% nii = make_nii(lfs,voxelSize);
	% save_nii(nii,'lfs_pdf.nii');

	% (2) RESHARP
	[lph,mask_ero] = resharp(unph,mask,voxelSize,ker_rad,tik_reg);
	nii = make_nii(mask_ero,voxelSize);
	save_nii(nii,'mask_ero.nii');
	nii = make_nii(lph,voxelSize);
	save_nii(nii,'lph.nii');


	% normalize to ppm unit
	TE = params.protocol_header.alTE{1}/1e6;
	B_0 = params.protocol_header.m_flMagneticFieldStrength;
	gamma = 2.675222e8;
	lfs = lph/(gamma*TE*B_0)*1e6; % unit ppm
	nii = make_nii(lfs,voxelSize);
	save_nii(nii,'lfs.nii');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% susceptibility inversion
	disp('--> (7/7) TV susceptibility inversion ...');

	% account for oblique slicing (head tilted)
	theta = -acos(params.protocol_header.sSliceArray.asSlice{1}.sNormal.dTra);

	tv_reg = 5e-4; % total variation regularization

	sus = tvdi(lfs,mask_ero,voxelSize,tv_reg,abs(img_cmb),theta);
	% sus_final = sus.*mask_ero;
	% nii = make_nii(sus_final,voxelSize);
	nii = make_nii(sus,voxelSize);
	save_nii(nii,'sus.nii');


	% for debugging purpose
	% save all the variables in all.mat
	save all.mat

end % end of the very first for loop

cd(init_dir);