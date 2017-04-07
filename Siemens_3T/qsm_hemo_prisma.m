function qsm_hemo_prisma(path_mag, path_ph, path_out, options)
%QSM_HEMO_PRISMA Quantitative susceptibility mapping from SWI sequence at PRISMA (3T) for hemorrhage.
%   QSM_HEMO_PRISMA(PATH_MAG, PATH_PH, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH_MAG     - directory of magnitude dicoms
%   PATH_PH      - directory of unfiltered phase dicoms
%   PATH_OUT     - directory to save nifti and/or matrixes   : QSM_SWI_PRISMA
%   OPTIONS      - parameter structure including fields below
%    .bet_thr    - threshold for BET brain mask              : 0.4
%    .bet_smooth - smoothness of BET brain mask at edges     : 2
%    .ph_unwrap  - 'prelude' or 'laplacian' or 'bestpath'    : 'bestpath'
%    .bkg_rm     - background field removal method(s)        : 'lbv'
%                  options: 'pdf','sharp','resharp','esharp','lbv'
%                  to try all e.g.: {'pdf','sharp','resharp','esharp','lbv'}
%    .t_svd      - truncation of SVD for SHARP               : 0.1
%    .smv_rad    - radius (mm) of SMV convolution kernel     : 3
%    .tik_reg    - Tikhonov regularization for resharp       : 1e-3
%    .cgs_num    - max interation number for RESHARP         : 500
%    .lbv_peel   - LBV layers to be peeled off               : 2
%    .lbv_tol    - LBV interation error tolerance            : 0.01
%    .tv_reg     - Total variation regularization parameter  : 5e-4
%    .tvdi_n     - iteration number of TVDI (nlcg)           : 500

if ~ exist('path_mag','var') || isempty(path_mag)
    error('Please input the directory of magnitude DICOMs')
end

if ~ exist('path_ph','var') || isempty(path_ph)
    error('Please input the directory of unfiltered phase DICOMs')
end

if ~ exist('path_out','var') || isempty(path_out)
    path_out = pwd;
    display('Current directory for output')
end

if ~ exist('options','var') || isempty(options)
    options = [];
end

if ~ isfield(options,'bet_thr')
    options.bet_thr = 0.4;
end

if ~ isfield(options,'bet_smooth')
    options.bet_smooth = 2;
end

if ~ isfield(options,'ph_unwrap')
    options.ph_unwrap = 'bestpath';
end

if ~ isfield(options,'bkg_rm')
    options.bkg_rm = 'lbv';
end

if ~ isfield(options,'t_svd')
    options.t_svd = 0.1;
end

if ~ isfield(options,'smv_rad')
    options.smv_rad = 3;
end

if ~ isfield(options,'tik_reg')
    options.tik_reg = 1e-3;
end

if ~ isfield(options,'cgs_num')
    options.cgs_num = 500;
end

if ~ isfield(options,'lbv_tol')
    options.lbv_tol = 0.01;
end

if ~ isfield(options,'lbv_peel')
    options.lbv_peel = 2;
end

if ~ isfield(options,'tv_reg')
    options.tv_reg = 5e-4;
end

if ~ isfield(options,'inv_num')
    options.inv_num = 500;
end


bet_thr    = options.bet_thr;
bet_smooth = options.bet_smooth;
ph_unwrap  = options.ph_unwrap;
bkg_rm     = options.bkg_rm;
t_svd      = options.t_svd;
smv_rad    = options.smv_rad;
tik_reg    = options.tik_reg;
cgs_num    = options.cgs_num;
lbv_tol    = options.lbv_tol;
lbv_peel   = options.lbv_peel;
tv_reg     = options.tv_reg;
inv_num    = options.inv_num;


% read in DICOMs of both magnitude and raw unfiltered phase images
% read in magnitude DICOMs
path_mag = cd(cd(path_mag));
mag_list = dir(path_mag);
mag_list = mag_list(~strncmpi('.', {mag_list.name}, 1));


for i = 1:length(mag_list)
    mag(:,:,i) = permute(single(dicomread([path_mag,filesep,mag_list(i).name])),[2,1]);
end

% get the sequence parameters
dicom_info = dicominfo([path_mag,filesep,mag_list(1).name]);
vox = [dicom_info.PixelSpacing(1), dicom_info.PixelSpacing(2), dicom_info.SliceThickness];
imsize = size(mag);

% angles!!! (z projections)
Xz = dicom_info.ImageOrientationPatient(3);
Yz = dicom_info.ImageOrientationPatient(6);
Zxyz = cross(dicom_info.ImageOrientationPatient(1:3),dicom_info.ImageOrientationPatient(4:6));
Zz = Zxyz(3);
%Zz = sqrt(1 - Xz^2 - Yz^2);
z_prjs = [Xz, Yz, Zz];

% read in phase DICOMs
path_ph = cd(cd(path_ph));
ph_list = dir(path_ph);
ph_list = ph_list(~strncmpi('.', {ph_list.name}, 1));

for i = 1:length(ph_list)
    ph(:,:,i) = permute(single(dicomread([path_ph,filesep,ph_list(i).name])),[2,1]);
    % covert to [-pi pi] range
    ph(:,:,i) = ph(:,:,i)/4095*2*pi - pi;
end


% define output directories
path_qsm = [path_out '/QSM_HEMO_PRISMA'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);


% save magnitude/phase in NIFTI form
mkdir('src');
nii = make_nii(mag,vox);
save_nii(nii,'src/mag.nii');
nii = make_nii(ph,vox);
save_nii(nii,'src/ph.nii');


% extract the brain and generate mask
setenv('bet_thr',num2str(bet_thr));
setenv('bet_smooth',num2str(bet_smooth));
[status,cmdout] = unix('rm BET*');
unix('bet2 src/mag.nii BET -f ${bet_thr} -m -w ${bet_smooth}');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);


% phase unwrapping, prelude is preferred!
if strcmpi('prelude',ph_unwrap)
    % unwrap phase with PRELUDE
    disp('--> unwrap aliasing phase ...');
    bash_script = ['prelude -a src/mag.nii -p src/ph.nii ' ...
        '-u unph.nii -m BET_mask.nii -n 8'];
    unix(bash_script);
    unix('gunzip -f unph.nii.gz');
    nii = load_nii('unph.nii');
    unph = double(nii.img);

elseif strcmpi('laplacian',ph_unwrap)
    % Ryan Topfer's Laplacian unwrapping
    Options.voxelSize = vox;
    unph = lapunwrap(ph, Options).*mask;
    nii = make_nii(unph, vox);
    save_nii(nii,'unph_lap.nii');

elseif strcmpi('bestpath',ph_unwrap)
    % unwrap the phase using best path
    [pathstr, ~, ~] = fileparts(which('3DSRNCP.m'));
    setenv('pathstr',pathstr);
    setenv('nv',num2str(imsize(1)));
    setenv('np',num2str(imsize(2)));
    setenv('ns',num2str(imsize(3)));

    fid = fopen('wrapped_phase.dat','w');
    fwrite(fid,ph,'float');
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
    unph = reshape(tmp - round(mean(tmp(mask==1))/(2*pi))*2*pi,imsize(1:3)).*mask;
    fclose(fid);

    nii = make_nii(unph,vox);
    save_nii(nii,'unph_bestpath.nii');

    fid = fopen('reliability.dat','r');
    reliability = fread(fid,'float');
    fclose(fid);

    % reliability_raw = fread(fid,'float');
    % reliability_raw = reshape(reliability_raw,imsize(1:3));
    reliability = reshape(reliability,imsize(1:3));
    reliability = 1./reliability.*mask;
    reliability(reliability <= 0.1) = 0;
    reliability(reliability > 0.1) = 1;

    nii = make_nii(reliability,vox);
    save_nii(nii,'reliability.nii');

else
    error('what unwrapping methods to use? prelude or laplacian or bestpath?')
end


% normalize total field (tfs) to ppm unit
tfs = unph/(2.675e8*dicom_info.EchoTime*dicom_info.MagneticFieldStrength)*1e9; % unit ppm

nii = make_nii(tfs,vox);
save_nii(nii,'tfs.nii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % background field removal
% % PDF
% if sum(strcmpi('pdf',bkg_rm))
%     disp('--> PDF to remove background field ...');
%     [lfs_pdf,mask_pdf] = projectionontodipolefields(tfs,mask,vox,smv_rad,mag,z_prjs);
%     % 3D 2nd order polyfit to remove any residual background
%     lfs_pdf= poly3d(lfs_pdf,mask_pdf);

%     % save nifti
%     mkdir('PDF');
%     nii = make_nii(lfs_pdf,vox);
%     save_nii(nii,'PDF/lfs_pdf.nii');

%     % inversion of susceptibility 
%     disp('--> TV susceptibility inversion on PDF...');
%     sus_pdf = tvdi(lfs_pdf,mask_pdf,vox,tv_reg,mag,z_prjs,inv_num); 

%     % save nifti
%     nii = make_nii(sus_pdf.*mask_pdf,vox);
%     save_nii(nii,'PDF/sus_pdf.nii');
% end

% % SHARP (t_svd: truncation threthold for t_svd)
% if sum(strcmpi('sharp',bkg_rm))
%     disp('--> SHARP to remove background field ...');
%     [lfs_sharp, mask_sharp] = sharp(tfs,mask,vox,smv_rad,t_svd);
%     % % 3D 2nd order polyfit to remove any residual background
%     % lfs_sharp= poly3d(lfs_sharp,mask_sharp);

%     % save nifti
%     mkdir('SHARP');
%     nii = make_nii(lfs_sharp,vox);
%     save_nii(nii,'SHARP/lfs_sharp.nii');
    
%     % inversion of susceptibility 
%     disp('--> TV susceptibility inversion on SHARP...');
%     sus_sharp = tvdi(lfs_sharp,mask_sharp,vox,tv_reg,mag,z_prjs,inv_num); 
   
%     % save nifti
%     nii = make_nii(sus_sharp.*mask_sharp,vox);
%     save_nii(nii,'SHARP/sus_sharp.nii');
% end

% % RE-SHARP (tik_reg: Tikhonov regularization parameter)
% if sum(strcmpi('resharp',bkg_rm))
%     disp('--> RESHARP to remove background field ...');
%     [lfs_resharp, mask_resharp] = resharp(tfs,mask,vox,smv_rad,tik_reg,cgs_num);
%     % % 3D 2nd order polyfit to remove any residual background
%     % lfs_resharp= poly3d(lfs_resharp,mask_resharp);

%     % save nifti
%     mkdir('RESHARP');
%     nii = make_nii(lfs_resharp,vox);
%     save_nii(nii,'RESHARP/lfs_resharp.nii');

%     % inversion of susceptibility 
%     disp('--> TV susceptibility inversion on RESHARP...');
%     sus_resharp = tvdi(lfs_resharp,mask_resharp,vox,tv_reg,mag,z_prjs,inv_num); 
   
%     % save nifti
%     nii = make_nii(sus_resharp.*mask_resharp,vox);
%     save_nii(nii,'RESHARP/sus_resharp.nii');
% end

% % E-SHARP (SHARP edge extension)
% if sum(strcmpi('esharp',bkg_rm))
%     disp('--> E-SHARP to remove background field ...');
%     Parameters.voxelSize             = vox; % in mm
%     Parameters.resharpRegularization = tik_reg ;
%     Parameters.resharpKernelRadius   = smv_rad ; % in mm
%     Parameters.radius                = [ 10 10 5 ] ;

%     % pad matrix size to even number
%     pad_size = mod(size(tfs),2);
%     tfs = padarray(tfs.*mask, pad_size, 'post');

%     % taking off additional 1 voxels from edge - not sure the outermost 
%     % phase data included in the original mask is reliable. 
%     mask_shaved = shaver( ( tfs ~= 0 ), 1 ) ; % 1 voxel taken off
%     totalField  = mask_shaved .* tfs ;

%     % resharp 
%     [reducedLocalField, maskReduced] = ...
%         resharp( totalField, ...
%                  double(mask_shaved), ...
%                  Parameters.voxelSize, ...
%                  Parameters.resharpKernelRadius, ...
%                  Parameters.resharpRegularization ) ;

%     % extrapolation ~ esharp 
%     reducedBackgroundField = maskReduced .* ( totalField - reducedLocalField) ;

%     extendedBackgroundField = extendharmonicfield( ...
%        reducedBackgroundField, mask, maskReduced, Parameters) ;

%     backgroundField = extendedBackgroundField + reducedBackgroundField ;
%     localField      = totalField - backgroundField ;

%     lfs_esharp      = localField(1+pad_size(1):end,1+pad_size(2):end,1+pad_size(3):end);
%     mask_esharp     = mask_shaved(1+pad_size(1):end,1+pad_size(2):end,1+pad_size(3):end);  

%     % % 3D 2nd order polyfit to remove any residual background
%     % lfs_esharp = poly3d(lfs_esharp,mask_esharp);

%     % save nifti
%     mkdir('ESHARP');
%     nii = make_nii(lfs_esharp,vox);
%     save_nii(nii,'ESHARP/lfs_esharp.nii');

%     % inversion of susceptibility 
%     disp('--> TV susceptibility inversion on ESHARP...');
%     sus_esharp = tvdi(lfs_esharp,mask_esharp,vox,tv_reg,mag,z_prjs,inv_num); 
   
%     % save nifti
%     nii = make_nii(sus_esharp.*mask_esharp,vox);
%     save_nii(nii,'ESHARP/sus_esharp.nii');
% end

% % LBV
% if sum(strcmpi('lbv',bkg_rm))
%     disp('--> LBV to remove background field ...');
%     lfs_lbv = LBV(tfs,mask,imsize,vox,lbv_tol,lbv_peel); % strip 2 layers
%     mask_lbv = ones(size(mask));
%     mask_lbv(lfs_lbv==0) = 0;
%     % 3D 2nd order polyfit to remove any residual background
%     lfs_lbv= poly3d(lfs_lbv,mask_lbv);

%     % save nifti
%     mkdir('LBV');
%     nii = make_nii(lfs_lbv,vox);
%     save_nii(nii,'LBV/lfs_lbv.nii');

%     % inversion of susceptibility 
%     disp('--> TV susceptibility inversion on lbv...');
%     sus_lbv = tvdi(lfs_lbv,mask_lbv,vox,tv_reg,mag,z_prjs,inv_num);   

%     % save nifti
%     nii = make_nii(sus_lbv.*mask_lbv,vox);
%     save_nii(nii,'LBV/sus_lbv.nii');
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lfs_lbv = LBV(tfs,mask,imsize,vox,lbv_tol,lbv_peel); % strip 2 layers
mask_lbv = ones(size(mask));
mask_lbv(lfs_lbv==0) = 0;
% 3D 2nd order polyfit to remove any residual background
lfs_lbv= lfs_lbv - poly3d(lfs_lbv,mask_lbv);

% save nifti
mkdir('LBV');
nii = make_nii(lfs_lbv,vox);
save_nii(nii,'LBV/lfs_lbv.nii');


% % (1) reliability weighted
% weights = abs(img_cmb)./max(abs(img_cmb(:))).*mask_lbv.*reliability;
% weights = smooth3(weights,'gaussian',[7,7,3],1);
% % weights = smooth3(weights,'gaussian',[7,7,3],0.5);
% nii = make_nii(weights,vox);
% save_nii(nii,'weights.nii');

% % inversion of susceptibility 
% disp('--> TV susceptibility inversion on lbv...');
% sus_lbv = tvdi(lfs_lbv,mask_lbv,vox,tv_reg, ...
%     weights,z_prjs,inv_num);   

% % save nifti
% nii = make_nii(sus_lbv.*mask_lbv,vox);
% save_nii(nii,'LBV/sus_lbv_weighted.nii');

% (2) magnitude threshold weighted
thr = 0.2;
img_smooth = smooth3(mag.*reliability,'gaussian',[7,7,3],2);
% calculate the median value of magnitude
tmp = img_smooth(logical(mask_lbv));
size(tmp)
median(tmp(:))
hemo_mask = mask_lbv;
hemo_mask((img_smooth<=thr*median(tmp(:)))) = 0;
nii = make_nii(hemo_mask,vox);
save_nii(nii,['hemo_mask_' num2str(thr) '.nii']);

weights = mag./max(mag(:)).*hemo_mask;

weights = smooth3(weights,'gaussian',[7,7,3],2);

nii = make_nii(weights,vox);
save_nii(nii,'weights.nii');



% inversion of susceptibility 
disp('--> TV susceptibility inversion on lbv...');
sus_lbv = tvdi(lfs_lbv,mask_lbv,vox,tv_reg, ...
    weights,z_prjs,inv_num);   

% save nifti
nii = make_nii(sus_lbv.*mask_lbv,vox);
save_nii(nii,'LBV/sus_lbv_masked.nii');



% super-position method to reduce artifacts
imsize = size(sus_lbv);
hemo_mask = ones(imsize);
% hemo_mask((weights<=0.1) & (mask_lbv>0)) = 0;
hemo_mask(sus_lbv>=0.5) = 0;

% % threshold the magnitude
% % untrusted worthy hemorrhage as 0, others 1
% hemo_mask = mask;
% hemo_mask(abs(img_smooth)<=0.3) = 0;

nii = make_nii(hemo_mask,vox);
save_nii(nii,'hemo_mask.nii');



% remove the hemorrhage dipole fields
% % LBV to remove hemorrhage
lfs_noHemo_lbv = LBV(lfs_lbv,mask_lbv.*hemo_mask,imsize,vox,0.0001,1); % strip 1 layers
% lfs_noHemo_lbv = LBV(lfs_lbv,1-(mask_lbv-hemo_mask_lbv),imsize,vox,0.01,1); % strip 1 layers
mask_noHemo_lbv = ones(imsize);
mask_noHemo_lbv(lfs_noHemo_lbv == 0) = 0;



% inversion on remaining part (noHemo) of the brain
weights = mag.*mask_lbv.*mask_noHemo_lbv;
sus_noHemo_lbv = tvdi(lfs_noHemo_lbv,mask_lbv.*mask_noHemo_lbv,vox,tv_reg,weights,z_prjs); 


% super add together
sus_super_lbv = sus_noHemo_lbv.*mask_lbv.*mask_noHemo_lbv + sus_lbv.*(mask_lbv - mask_noHemo_lbv);


% baseline fix
Nx = size(sus_super_lbv,1);
Ny = size(sus_super_lbv,2);
Nz = size(sus_super_lbv,3);
FOV = vox.*[Nx,Ny,Nz];
FOVx = FOV(1);
FOVy = FOV(2);
FOVz = FOV(3);

x = -Nx/2:Nx/2-1;
y = -Ny/2:Ny/2-1;
z = -Nz/2:Nz/2-1;
[kx,ky,kz] = ndgrid(x/FOVx,y/FOVy,z/FOVz);
D = 1/3 - (kx.*z_prjs(1)+ky.*z_prjs(2)+kz.*z_prjs(3)).^2./(kx.^2 + ky.^2 + kz.^2);
D(floor(Nx/2+1),floor(Ny/2+1),floor(Nz/2+1)) = 0;
D = fftshift(D);

W = weights;

x1 = W.*ifftn(D.*fftn(mask_lbv - mask_noHemo_lbv));
x2 = W.*(lfs_lbv - ifftn(D.*fftn(sus_super_lbv)));
x1 = x1(:);
x2 = x2(:);
o = real(x1'*x2/(x1'*x1))

sus_super_fix_lbv = sus_super_lbv + o*(mask_lbv - mask_noHemo_lbv);
nii = make_nii(sus_super_fix_lbv,vox);
save_nii(nii,'LBV/sus_super_fix_lbv.nii');



save('all.mat','-v7.3');
cd(init_dir);



