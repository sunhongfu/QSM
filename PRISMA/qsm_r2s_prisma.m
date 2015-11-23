function qsm_r2s_prisma(path_mag, path_ph, path_out, options)
%QSM_R2S_PRISMA Quantitative susceptibility mapping from R2s sequence at PRISMA (3T).
%   QSM_R2S_PRISMA(PATH_MAG, PATH_PH, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH_MAG    - directory of magnitude dicoms
%   PATH_PH     - directory of unfiltered phase dicoms
%   PATH_OUT    - directory to save nifti and/or matrixes   : QSM_SWI_PRISMA
%   OPTIONS     - parameter structure including fields below
%    .readout   - multi-echo 'unipolar' or 'bipolar'        : 'unipolar'
%    .r_mask    - whether to enable the extra masking       : 1
%    .fit_thr   - extra filtering based on the fit residual : 40
%    .bet_thr   - threshold for BET brain mask              : 0.4
%    .ph_unwrap - 'prelude' or 'bestpath'                   : 'prelude'
%    .bkg_rm    - background field removal method(s)        : 'resharp'
%                 options: 'pdf','sharp','resharp','esharp','lbv'
%                 to try all e.g.: {'pdf','sharp','resharp','esharp','lbv'}
%    .t_svd     - truncation of SVD for SHARP               : 0.1
%    .smv_rad   - radius (mm) of SMV convolution kernel     : 3
%    .tik_reg   - Tikhonov regularization for resharp       : 1e-3
%    .cgs_num   - max interation number for RESHARP         : 500
%    .lbv_peel  - LBV layers to be peeled off               : 2
%    .lbv_tol   - LBV interation error tolerance            : 0.01
%    .tv_reg    - Total variation regularization parameter  : 5e-4
%    .tvdi_n    - iteration number of TVDI (nlcg)           : 500

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

if ~ isfield(options,'readout')
    options.readout = 'unipolar';
end

if ~ isfield(options,'r_mask')
    options.r_mask = 1;
end

if ~ isfield(options,'fit_thr')
    options.fit_thr = 40;
end

if ~ isfield(options,'bet_thr')
    options.bet_thr = 0.4;
end

if ~ isfield(options,'ph_unwrap')
    options.ph_unwrap = 'prelude';
end

if ~ isfield(options,'bkg_rm')
    options.bkg_rm = 'resharp';
    % options.bkg_rm = {'pdf','sharp','resharp','esharp','lbv'};
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



readout   = options.readout;
r_mask    = options.r_mask;
fit_thr   = options.fit_thr;
bet_thr   = options.bet_thr;
ph_unwrap = options.ph_unwrap;
bkg_rm    = options.bkg_rm;
t_svd     = options.t_svd;
smv_rad   = options.smv_rad;
tik_reg   = options.tik_reg;
cgs_num   = options.cgs_num;
lbv_tol   = options.lbv_tol;
lbv_peel  = options.lbv_peel;
tv_reg    = options.tv_reg;
inv_num   = options.inv_num;


% read in DICOMs of both magnitude and raw unfiltered phase images
% read in magnitude DICOMs
path_mag = cd(cd(path_mag));
mag_list = dir(path_mag);

% get the sequence parameters
dicom_info = dicominfo([path_mag,filesep,mag_list(3).name]);
EchoTrainLength = dicom_info.EchoTrainLength;
for i = 1:EchoTrainLength % read in TEs
    dicom_info = dicominfo([path_mag,filesep,mag_list(3+(i-1)*(length(mag_list)-2)./EchoTrainLength).name]);
    TE(dicom_info.EchoNumber) = dicom_info.EchoTime*1e-3;
end
vox = [dicom_info.PixelSpacing(1), dicom_info.PixelSpacing(2), dicom_info.SliceThickness];

% angles!!! (z projections)
Xz = dicom_info.ImageOrientationPatient(3);
Yz = dicom_info.ImageOrientationPatient(6);
Zz = sqrt(1 - Xz^2 - Yz^2);
z_prjs = [Xz, Yz, Zz];

for i = 3:length(mag_list)
    [NS,NE] = ind2sub([(length(mag_list)-2)./EchoTrainLength,EchoTrainLength],i-2);
    mag(:,:,NS,NE) = permute(single(dicomread([path_mag,filesep,mag_list(i).name])),[2,1]);
end
% size of matrix
imsize = size(mag);

% read in phase DICOMs
path_ph = cd(cd(path_ph));
ph_list = dir(path_ph);

for i = 3:length(ph_list)
    [NS,NE] = ind2sub([(length(ph_list)-2)./EchoTrainLength,EchoTrainLength],i-2);
    ph(:,:,NS,NE) = permute(single(dicomread([path_ph,filesep,ph_list(i).name])),[2,1]);    % covert to [-pi pi] range
    ph(:,:,NS,NE) = ph(:,:,NS,NE)/4095*2*pi - pi;
end


% define output directories
path_qsm = [path_out '/QSM_R2S_PRISMA'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);


% save magnitude and raw phase niftis for each echo
mkdir('src')
for echo = 1:imsize(4)
    nii = make_nii(mag(:,:,:,echo),vox);
    save_nii(nii,['src/mag' num2str(echo) '.nii']);
    nii = make_nii(ph(:,:,:,echo),vox);
    save_nii(nii,['src/ph' num2str(echo) '.nii']);
end


% brain extraction
% generate mask from magnitude of the 1th echo
disp('--> extract brain volume and generate mask ...');
setenv('bet_thr',num2str(bet_thr));
[status,cmdout] = unix('rm BET*');
unix('bet src/mag1.nii BET -f ${bet_thr} -m -R');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);


% phase offset correction
% if unipolar
if strcmpi('unipolar',readout)
    ph_corr = geme_cmb(mag.*exp(1j*ph),vox,TE,mask);
% if bipolar
elseif strcmpi('bipolar',readout)
    ph_corr = zeros(imsize);
    ph_corr(:,:,:,1:2:end) = geme_cmb(mag(:,:,:,1:2:end).*exp(1j*ph(:,:,:,1:2:end)),vox,TE(1:2:end),mask);
    ph_corr(:,:,:,2:2:end) = geme_cmb(mag(:,:,:,2:2:end).*exp(1j*ph(:,:,:,2:2:end)),vox,TE(2:2:end),mask);
else
    error('is the sequence unipolar or bipolar readout?')
end

% save offset corrected phase niftis
for echo = 1:imsize(4)
    nii = make_nii(ph_corr(:,:,:,echo),vox);
    save_nii(nii,['src/ph_corr' num2str(echo) '.nii']);
end


% unwrap phase from each echo
if strcmpi('prelude',ph_unwrap)
    disp('--> unwrap aliasing phase for all TEs using prelude...');
    setenv('echo_num',num2str(imsize(4)));
    bash_command = sprintf(['for ph in src/ph_corr[1-$echo_num].nii\n' ...
    'do\n' ...
    '	base=`basename $ph`;\n' ...
    '	dir=`dirname $ph`;\n' ...
    '	mag=$dir/"mag"${base:7};\n' ...
    '	unph="unph"${base:7};\n' ...
    '	prelude -a $mag -p $ph -u $unph -m BET_mask.nii -n 12&\n' ...
    'done\n' ...
    'wait\n' ...
    'gunzip -f unph*.gz\n']);
    unix(bash_command);

    unph = zeros(imsize);
    for echo = 1:imsize(4)
        nii = load_nii(['unph' num2str(echo) '.nii']);
        unph(:,:,:,echo) = double(nii.img);
    end


elseif strcmpi('bestpath',ph_unwrap)
    % unwrap the phase using best path
    disp('--> unwrap aliasing phase using bestpath...');
    mask_unwrp = uint8(abs(mask)*255);
    fid = fopen('mask_unwrp.dat','w');
    fwrite(fid,mask_unwrp,'uchar');
    fclose(fid);

    [pathstr, ~, ~] = fileparts(which('3DSRNCP.m'));
    setenv('pathstr',pathstr);
    setenv('nv',num2str(imsize(1)));
    setenv('np',num2str(imsize(2)));
    setenv('ns',num2str(imsize(3)));

    unph = zeros(imsize);

    for echo_num = 1:imsize(4)
        setenv('echo_num',num2str(echo_num));
        fid = fopen(['wrapped_phase' num2str(echo_num) '.dat'],'w');
        fwrite(fid,ph_corr(:,:,:,echo_num),'float');
        fclose(fid);

        bash_script = ['${pathstr}/3DSRNCP wrapped_phase${echo_num}.dat mask_unwrp.dat ' ...
            'unwrapped_phase${echo_num}.dat $nv $np $ns reliability${echo_num}.dat'];
        unix(bash_script) ;

        fid = fopen(['unwrapped_phase' num2str(echo_num) '.dat'],'r');
        tmp = fread(fid,'float');
        % tmp = tmp - tmp(1);
        unph(:,:,:,echo_num) = reshape(tmp - round(mean(tmp(mask==1))/(2*pi))*2*pi ,imsize(1:3)).*mask;
        fclose(fid);

        fid = fopen(['reliability' num2str(echo_num) '.dat'],'r');
        reliability_raw = fread(fid,'float');
        reliability_raw = reshape(reliability_raw,imsize(1:3));
        fclose(fid);

        nii = make_nii(reliability_raw.*mask,vox);
        save_nii(nii,['reliability_raw' num2str(echo_num) '.nii']);
    end

    nii = make_nii(unph,vox);
    save_nii(nii,'unph_bestpath.nii');

else
    error('what unwrapping methods to use? prelude or bestpath?')
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

for echo = 2:imsize(4)
    meandiff = unph(:,:,:,echo)-unph(:,:,:,1)-double(echo-1)*unph_diff;
    meandiff = meandiff(mask==1);
    meandiff = mean(meandiff(:));
    njump = round(meandiff/(2*pi));
    disp(['    ' num2str(njump) ' 2pi jumps for TE' num2str(echo)]);
    unph(:,:,:,echo) = unph(:,:,:,echo) - njump*2*pi;
    unph(:,:,:,echo) = unph(:,:,:,echo).*mask;
end


% fit phase images with echo times
disp('--> magnitude weighted LS fit of phase to TE ...');
[tfs, fit_residual] = echofit(unph,mag,TE,0); 


% extra filtering according to fitting residuals
if r_mask
    % generate reliability map
    fit_residual_blur = smooth3(fit_residual,'box',round(1./vox)*2+1); 
    nii = make_nii(fit_residual_blur,vox);
    save_nii(nii,'fit_residual_blur.nii');
    R = ones(size(fit_residual_blur));
    R(fit_residual_blur >= fit_thr) = 0;
else
    R = 1;
end


% normalize to main field
% ph = gamma*dB*TE
% dB/B = ph/(gamma*TE*B0)
% units: TE s, gamma 2.675e8 rad/(sT), B0 4.7T
tfs = tfs/(2.675e8*dicom_info.MagneticFieldStrength)*1e6; % unit ppm

nii = make_nii(tfs,vox);
save_nii(nii,'tfs.nii');


% background field removal
% PDF
if sum(strcmpi('pdf',bkg_rm))
    disp('--> PDF to remove background field ...');
    [lfs_pdf,mask_pdf] = projectionontodipolefields(tfs,mask.*R,vox,smv_rad,mag(:,:,:,end),z_prjs);
    % 3D 2nd order polyfit to remove any residual background
    lfs_pdf= poly3d(lfs_pdf,mask_pdf);

    % save nifti
    mkdir('PDF');
    nii = make_nii(lfs_pdf,vox);
    save_nii(nii,'PDF/lfs_pdf.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on PDF...');
    sus_pdf = tvdi(lfs_pdf,mask_pdf,vox,tv_reg,mag(:,:,:,end),z_prjs,inv_num); 

    % save nifti
    nii = make_nii(sus_pdf.*mask_pdf,vox);
    save_nii(nii,'PDF/sus_pdf.nii');
end

% SHARP (t_svd: truncation threthold for t_svd)
if sum(strcmpi('sharp',bkg_rm))
    disp('--> SHARP to remove background field ...');
    [lfs_sharp, mask_sharp] = sharp(tfs,mask.*R,vox,smv_rad,t_svd);
    % % 3D 2nd order polyfit to remove any residual background
    % lfs_sharp= poly3d(lfs_sharp,mask_sharp);

    % save nifti
    mkdir('SHARP');
    nii = make_nii(lfs_sharp,vox);
    save_nii(nii,'SHARP/lfs_sharp.nii');
    
    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on SHARP...');
    sus_sharp = tvdi(lfs_sharp,mask_sharp,vox,tv_reg,mag(:,:,:,end),z_prjs,inv_num); 
   
    % save nifti
    nii = make_nii(sus_sharp.*mask_sharp,vox);
    save_nii(nii,'SHARP/sus_sharp.nii');
end

% RE-SHARP (tik_reg: Tikhonov regularization parameter)
if sum(strcmpi('resharp',bkg_rm))
    disp('--> RESHARP to remove background field ...');
    [lfs_resharp, mask_resharp] = resharp(tfs,mask.*R,vox,smv_rad,tik_reg,cgs_num);
    % % 3D 2nd order polyfit to remove any residual background
    % lfs_resharp= poly3d(lfs_resharp,mask_resharp);

    % save nifti
    mkdir('RESHARP');
    nii = make_nii(lfs_resharp,vox);
    save_nii(nii,'RESHARP/lfs_resharp.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on RESHARP...');
    sus_resharp = tvdi(lfs_resharp,mask_resharp,vox,tv_reg,mag(:,:,:,end),z_prjs,inv_num); 
   
    % save nifti
    nii = make_nii(sus_resharp.*mask_resharp,vox);
    save_nii(nii,'RESHARP/sus_resharp.nii');
end

% E-SHARP (SHARP edge extension)
if sum(strcmpi('esharp',bkg_rm))
    disp('--> E-SHARP to remove background field ...');
    Parameters.voxelSize             = vox; % in mm
    Parameters.resharpRegularization = tik_reg ;
    Parameters.resharpKernelRadius   = smv_rad ; % in mm
    Parameters.radius                = [ 10 10 5 ] ;

    % pad matrix size to even number
    pad_size = mod(size(tfs),2);
    tfs = tfs.*mask.*R;
    tfs = padarray(tfs, pad_size, 'post');

    % taking off additional 1 voxels from edge - not sure the outermost 
    % phase data included in the original mask is reliable. 
    mask_shaved = shaver( ( tfs ~= 0 ), 1 ) ; % 1 voxel taken off
    totalField  = mask_shaved .* tfs ;

    % resharp 
    [reducedLocalField, maskReduced] = ...
        resharp( totalField, ...
                 double(mask_shaved), ...
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
    mask_esharp     = mask_shaved(1+pad_size(1):end,1+pad_size(2):end,1+pad_size(3):end);  

    % % 3D 2nd order polyfit to remove any residual background
    % lfs_esharp = poly3d(lfs_esharp,mask_esharp);

    % save nifti
    mkdir('ESHARP');
    nii = make_nii(lfs_esharp,vox);
    save_nii(nii,'ESHARP/lfs_esharp.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on ESHARP...');
    sus_esharp = tvdi(lfs_esharp,mask_esharp,vox,tv_reg,mag(:,:,:,end),z_prjs,inv_num); 
   
    % save nifti
    nii = make_nii(sus_esharp.*mask_esharp,vox);
    save_nii(nii,'ESHARP/sus_esharp.nii');
end

% LBV
if sum(strcmpi('lbv',bkg_rm))
    disp('--> LBV to remove background field ...');
    lfs_lbv = LBV(tfs,mask.*R,imsize(1:3),vox,lbv_tol,lbv_peel); % strip 2 layers
    mask_lbv = ones(imsize(1:3));
    mask_lbv(lfs_lbv==0) = 0;
    % 3D 2nd order polyfit to remove any residual background
    lfs_lbv= poly3d(lfs_lbv,mask_lbv);

    % save nifti
    mkdir('LBV');
    nii = make_nii(lfs_lbv,vox);
    save_nii(nii,'LBV/lfs_lbv.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on lbv...');
    sus_lbv = tvdi(lfs_lbv,mask_lbv,vox,tv_reg,mag(:,:,:,end),z_prjs,inv_num);   

    % save nifti
    nii = make_nii(sus_lbv.*mask_lbv,vox);
    save_nii(nii,'LBV/sus_lbv.nii');
end


save('all.mat','-v7.3');
cd(init_dir);

