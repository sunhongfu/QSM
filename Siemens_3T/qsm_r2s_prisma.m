function qsm_r2s_prisma(path_mag, path_ph, path_out, options)
%QSM_R2S_PRISMA Quantitative susceptibility mapping from R2s sequence at PRISMA (3T).
%   QSM_R2S_PRISMA(PATH_MAG, PATH_PH, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH_MAG     - directory of magnitude dicoms
%   PATH_PH      - directory of unfiltered phase dicoms
%   PATH_OUT     - directory to save nifti and/or matrixes   : QSM_SWI_PRISMA
%   OPTIONS      - parameter structure including fields below
%    .readout    - multi-echo 'unipolar' or 'bipolar'        : 'unipolar'
%    .r_mask     - whether to enable the extra masking       : 1
%    .fit_thr    - extra filtering based on the fit residual : 40
%    .bet_thr    - threshold for BET brain mask              : 0.4
%    .bet_smooth - smoothness of BET brain mask at edges     : 2
%    .ph_unwrap  - 'prelude' or 'bestpath'                   : 'prelude'
%    .bkg_rm     - background field removal method(s)        : 'resharp'
%                  options: 'pdf','sharp','resharp','esharp','lbv'
%                  to try all e.g.: {'pdf','sharp','resharp','esharp','lbv'}
%                  set to 'lnqsm' or 'tfi' for LN-QSM and TFI methods
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
    disp('Current directory for output')
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

if ~ isfield(options,'bet_smooth')
    options.bet_smooth = 2;
end

if ~ isfield(options,'ph_unwrap')
    options.ph_unwrap = 'bestpath';
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

if ~ isfield(options,'interp')
    options.interp = 0;
end

readout    = options.readout;
r_mask     = options.r_mask;
fit_thr    = options.fit_thr;
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
interp     = options.interp;

% read in DICOMs of both magnitude and raw unfiltered phase images
% read in magnitude DICOMs
path_mag = cd(cd(path_mag));
mag_list = dir(path_mag);
mag_list = mag_list(~strncmpi('.', {mag_list.name}, 1));

% get the sequence parameters
dicom_info = dicominfo([path_mag,filesep,mag_list(end).name]);
% EchoTrainLength = dicom_info.EchoTrainLength;
EchoTrainLength = dicom_info.EchoNumbers;
for i = 1:EchoTrainLength % read in TEs
    dicom_info = dicominfo([path_mag,filesep,mag_list(1+(i-1)*(length(mag_list))./EchoTrainLength).name]);
    TE(dicom_info.EchoNumbers) = dicom_info.EchoTime*1e-3;
end
vox = [dicom_info.PixelSpacing(1), dicom_info.PixelSpacing(2), dicom_info.SliceThickness];

% angles!!! (z projections)
Xz = dicom_info.ImageOrientationPatient(3);
Yz = dicom_info.ImageOrientationPatient(6);
%Zz = sqrt(1 - Xz^2 - Yz^2);
Zxyz = cross(dicom_info.ImageOrientationPatient(1:3),dicom_info.ImageOrientationPatient(4:6));
Zz = Zxyz(3);
z_prjs = [Xz, Yz, Zz];

for i = 1:length(mag_list)
    [NS,NE] = ind2sub([length(mag_list)./EchoTrainLength,EchoTrainLength],i);
    mag(:,:,NS,NE) = permute(single(dicomread([path_mag,filesep,mag_list(i).name])),[2,1]);
end
% size of matrix
imsize = size(mag);

% read in phase DICOMs
path_ph = cd(cd(path_ph));
ph_list = dir(path_ph);
ph_list = ph_list(~strncmpi('.', {ph_list.name}, 1));

for i = 1:length(ph_list)
    [NS,NE] = ind2sub([length(ph_list)./EchoTrainLength,EchoTrainLength],i);
    ph(:,:,NS,NE) = permute(single(dicomread([path_ph,filesep,ph_list(i).name])),[2,1]);    % covert to [-pi pi] range
    ph(:,:,NS,NE) = ph(:,:,NS,NE)/4095*2*pi - pi;
end



% interpolation into isotropic
if interp
    img = mag.*exp(1j*ph);
    k = fftshift(fftshift(fftshift(fft(fft(fft(img,[],1),[],2),[],3),1),2),3);
    % find the finest resolution
    minvox = min(vox);
    % update matrix size
    pad_size =  round((vox.*imsize(1:3)/minvox - imsize(1:3))/2);
    k = padarray(k, pad_size);
    img = ifft(ifft(ifft(ifftshift(ifftshift(ifftshift(k,1),2),3),[],1),[],2),[],3);
    clear k;
    imsize_old = imsize;
    imsize = size(img);
    vox = imsize_old(1:3).*vox./imsize(1:3);
    mag = abs(img);
    ph = angle(img);
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
setenv('bet_smooth',num2str(bet_smooth));
[~,~] = unix('rm BET*');
unix('bet2 src/mag1.nii BET -f ${bet_thr} -m -w ${bet_smooth}');
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


% LN-QSM
if sum(strcmpi('lnqsm',bkg_rm))
    mkdir LN-QSM

    maskR = mask.*R;
    iMag = sqrt(sum(mag.^2,4));
    Res_wt = iMag.*maskR;
    Res_wt = Res_wt/sum(Res_wt(:))*sum(maskR(:));

    % mask_erosion
    r = 1; 
    [X,Y,Z] = ndgrid(-r:r,-r:r,-r:r);
    h = (X.^2/r^2 + Y.^2/r^2 + Z.^2/r^2 <= 1);
    ker = h/sum(h(:));
    imsize = size(mask);
    mask_tmp = convn(mask,ker,'same');
    mask_ero1 = zeros(imsize);
    mask_ero1(mask_tmp > 0.999999) = 1; % no error tolerance
    r = 2; 
    [X,Y,Z] = ndgrid(-r:r,-r:r,-r:r);
    h = (X.^2/r^2 + Y.^2/r^2 + Z.^2/r^2 <= 1);
    ker = h/sum(h(:));
    imsize = size(mask);
    mask_tmp = convn(mask,ker,'same');
    mask_ero2 = zeros(imsize);
    mask_ero2(mask_tmp > 0.999999) = 1; % no error tolerance
    r = 3; 
    [X,Y,Z] = ndgrid(-r:r,-r:r,-r:r);
    h = (X.^2/r^2 + Y.^2/r^2 + Z.^2/r^2 <= 1);
    ker = h/sum(h(:));
    imsize = size(mask);
    mask_tmp = convn(mask,ker,'same');
    mask_ero3 = zeros(imsize);
    mask_ero3(mask_tmp > 0.999999) = 1; % no error tolerance
    
    % (1) no erosion, keep the full brain
    P = maskR + 30*(1 - maskR);
    chi_ero0_500 = tikhonov_qsm(tfs, Res_wt.*maskR, 1, maskR, maskR, 0, 1e-4, 0.001, 0, vox, P, z_prjs, 500);
    nii = make_nii(chi_ero0_500.*maskR,vox);
    save_nii(nii,'LN-QSM/chi_ero0_tik_1e-3_tv_1e-4_500.nii');

    % (2) erode 1 voxel from brain edge
    P = mask_ero1 + 30*(1 - mask_ero1);
    chi_ero1_500 = tikhonov_qsm(tfs, Res_wt.*mask_ero1, 1, mask_ero1, mask_ero1, 0, 1e-4, 0.001, 0, vox, P, z_prjs, 500);
    nii = make_nii(chi_ero1_500.*mask_ero1,vox);
    save_nii(nii,'LN-QSM/chi_ero1_tik_1e-3_tv_1e-4_500.nii');

    % (3) erode 2 voxel from brain edge
    P = mask_ero2 + 30*(1 - mask_ero2);
    chi_ero2_500 = tikhonov_qsm(tfs, Res_wt.*mask_ero2, 1, mask_ero2, mask_ero2, 0, 1e-4, 0.001, 0, vox, P, z_prjs, 500);
    nii = make_nii(chi_ero2_500.*mask_ero2,vox);
    save_nii(nii,'LN-QSM/chi_ero2_tik_1e-3_tv_1e-4_500.nii');

    % (4) erode 3 voxel from brain edge
    P = mask_ero3 + 30*(1 - mask_ero3);
    chi_ero3_500 = tikhonov_qsm(tfs, Res_wt.*mask_ero3, 1, mask_ero3, mask_ero3, 0, 1e-4, 0.001, 0, vox, P, z_prjs, 500);
    nii = make_nii(chi_ero3_500.*mask_ero3,vox);
    save_nii(nii,'LN-QSM/chi_ero3_tik_1e-3_tv_1e-4_500.nii');
    
%     % pad zeros
%     tfs_pad = padarray(tfs,[0 0 80]);
%     Res_wt_pad = padarray(Res_wt,[0 0 80]);
%     
%     mask_ero0_pad = padarray(maskR,[0 0 80]);
%     P_pad = mask_ero0_pad + 30*(1 - mask_ero0_pad);
%     chi_ero0_500_pad = tikhonov_qsm(tfs_pad, Res_wt_pad.*mask_ero0_pad, 1, mask_ero0_pad, mask_ero0_pad, 0, 1e-4, 0.001, 0, vox, P_pad, z_prjs, 500);
%     nii = make_nii(chi_ero0_500_pad(:,:,81:end-80).*maskR,vox);
%     save_nii(nii,'LN-QSM/chi_ero0_tik_1e-3_tv_1e-4_500_pad80.nii');
%     
%     mask_ero3_pad = padarray(mask_ero3,[0 0 80]);
%     P_pad = mask_ero3_pad + 30*(1 - mask_ero3_pad);
%     chi_ero3_500_pad = tikhonov_qsm(tfs_pad, Res_wt_pad.*mask_ero3_pad, 1, mask_ero3_pad, mask_ero3_pad, 0, 1e-4, 0.001, 0, vox, P_pad, z_prjs, 500);
%     nii = make_nii(chi_ero3_500_pad(:,:,81:end-80).*mask_ero3,vox);
%     save_nii(nii,'LN-QSM/chi_ero3_tik_1e-3_tv_1e-4_500_pad80.nii');
end

    
% background field removal
% PDF
if sum(strcmpi('pdf',bkg_rm))
    disp('--> PDF to remove background field ...');
    [lfs_pdf,mask_pdf] = projectionontodipolefields(tfs,mask.*R,vox,smv_rad,mag(:,:,:,end),z_prjs);
    % 3D 2nd order polyfit to remove any residual background
    lfs_pdf= lfs_pdf - poly3d(lfs_pdf,mask_pdf);

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
    % iLSQR
    chi_iLSQR = QSM_iLSQR(lfs_resharp*(2.675e8*dicom_info.MagneticFieldStrength)/1e6,mask_resharp,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
    nii = make_nii(chi_iLSQR,vox);
    save_nii(nii,['RESHARP/chi_iLSQR_smvrad' num2str(smv_rad) '.nii']);
    
    % % MEDI
    % %%%%% normalize signal intensity by noise to get SNR %%%
    % %%%% Generate the Magnitude image %%%%
    % iMag = sqrt(sum(mag.^2,4));
    % % [iFreq_raw N_std] = Fit_ppm_complex(ph_corr);
    % matrix_size = single(imsize(1:3));
    % voxel_size = vox;
    % delta_TE = TE(2) - TE(1);
    % B0_dir = z_prjs';
    % CF = dicom_info.ImagingFrequency *1e6;
    % iFreq = [];
    % N_std = 1;
    % RDF = lfs_resharp*2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6;
    % Mask = mask_resharp;
    % save RDF.mat RDF iFreq iMag N_std Mask matrix_size...
    %      voxel_size delta_TE CF B0_dir;
    % QSM = MEDI_L1('lambda',1000);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,['RESHARP/MEDI1000_RESHARP_smvrad' num2str(smv_rad) '.nii']);
    % QSM = MEDI_L1('lambda',2000);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,['RESHARP/MEDI2000_RESHARP_smvrad' num2str(smv_rad) '.nii']);
    % QSM = MEDI_L1('lambda',1500);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,['RESHARP/MEDI1500_RESHARP_smvrad' num2str(smv_rad) '.nii']);
    % QSM = MEDI_L1('lambda',5000);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,['RESHARP/MEDI5000_RESHARP_smvrad' num2str(smv_rad) '.nii']);
    % QSM = MEDI_L1('lambda',500);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,['RESHARP/MEDI500_RESHARP_smvrad' num2str(smv_rad) '.nii']);

    % TVDI method
    sus_resharp = tvdi(lfs_resharp,mask_resharp,vox,tv_reg,mag(:,:,:,end),z_prjs,inv_num); 
    nii = make_nii(sus_resharp.*mask_resharp,vox);
    save_nii(nii,'RESHARP/sus_resharp.nii');

end

% V-SHARP
if sum(strcmpi('vsharp',bkg_rm))
    disp('--> V-SHARP to remove background field ...');
    smvsize = 12;
    [lfs_vsharp, mask_vsharp] = V_SHARP(tfs ,single(mask.*R),'smvsize',smvsize,'voxelsize',vox);
    
    % save nifti
    mkdir('VSHARP');
    nii = make_nii(lfs_vsharp,vox);
    save_nii(nii,'VSHARP/lfs_vsharp.nii');
    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on RESHARP...');
    sus_vsharp = tvdi(lfs_vsharp,mask_vsharp,vox,tv_reg,mag(:,:,:,end),z_prjs,inv_num); 
   
    % save nifti
    nii = make_nii(sus_vsharp.*mask_vsharp,vox);
    save_nii(nii,'VSHARP/sus_vsharp.nii');
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
    lfs_lbv = lfs_lbv - poly3d(lfs_lbv,mask_lbv);

    % save nifti
    mkdir('LBV');
    nii = make_nii(lfs_lbv,vox);
    save_nii(nii,'LBV/lfs_lbv.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on LBV...');
    % iLSQR
    chi_iLSQR = QSM_iLSQR(lfs_lbv*(2.675e8*dicom_info.MagneticFieldStrength)/1e6,mask_lbv,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
    nii = make_nii(chi_iLSQR,vox);
    save_nii(nii,['LBV/chi_iLSQR_smvrad' num2str(smv_rad) '.nii']);
    
    % MEDI
    %%%%% normalize signal intensity by noise to get SNR %%%
    %%%% Generate the Magnitude image %%%%
    iMag = sqrt(sum(mag.^2,4));
    % [iFreq_raw N_std] = Fit_ppm_complex(ph_corr);
    matrix_size = single(imsize(1:3));
    voxel_size = vox;
    delta_TE = TE(2) - TE(1);
    B0_dir = z_prjs';
    CF = dicom_info.ImagingFrequency *1e6;
    iFreq = [];
    N_std = 1;
    RDF = lfs_lbv*2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6;
    Mask = mask_lbv;
    save RDF.mat RDF iFreq iMag N_std Mask matrix_size...
            voxel_size delta_TE CF B0_dir;
    QSM = MEDI_L1('lambda',1000);
    nii = make_nii(QSM.*Mask,vox);
    save_nii(nii,['LBV/MEDI1000_lbv_smvrad' num2str(smv_rad) '.nii']);
    QSM = MEDI_L1('lambda',2000);
    nii = make_nii(QSM.*Mask,vox);
    save_nii(nii,['LBV/MEDI2000_lbv_smvrad' num2str(smv_rad) '.nii']);
    QSM = MEDI_L1('lambda',1500);
    nii = make_nii(QSM.*Mask,vox);
    save_nii(nii,['LBV/MEDI1500_lbv_smvrad' num2str(smv_rad) '.nii']);
    QSM = MEDI_L1('lambda',5000);
    nii = make_nii(QSM.*Mask,vox);
    save_nii(nii,['LBV/MEDI5000_lbv_smvrad' num2str(smv_rad) '.nii']);
    QSM = MEDI_L1('lambda',500);
    nii = make_nii(QSM.*Mask,vox);
    save_nii(nii,['LBV/MEDI500_lbv_smvrad' num2str(smv_rad) '.nii']);

    % TVDI method
    sus_lbv = tvdi(lfs_lbv,mask_lbv,vox,tv_reg,mag(:,:,:,end),z_prjs,inv_num); 
    nii = make_nii(sus_lbv.*mask_lbv,vox);
    save_nii(nii,'LBV/sus_lbv.nii');
end



% TFI
if sum(strcmpi('tfi',bkg_rm))
    mkdir TFI
    %%%% Generate the Magnitude image %%%%
    iMag = sqrt(sum(mag.^2,4));
    matrix_size = single(imsize(1:3));
    voxel_size = vox;
    delta_TE = TE(2) - TE(1);
    B0_dir = z_prjs';
    CF = dicom_info.ImagingFrequency *1e6;
    N_std = 1;

    % (4) TFI of 3 voxels erosion
    iFreq = tfs*2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6;
    % erode the mask (full mask to 3mm erosion)
    % apply R
    mask = mask.*R;
    % mask_erosion
    r = 1; 
    [X,Y,Z] = ndgrid(-r:r,-r:r,-r:r);
    h = (X.^2/r^2 + Y.^2/r^2 + Z.^2/r^2 <= 1);
    ker = h/sum(h(:));
    imsize = size(mask);
    mask_tmp = convn(mask,ker,'same');
    mask_ero1 = zeros(imsize);
    mask_ero1(mask_tmp > 0.999999) = 1; % no error tolerance
    r = 2; 
    [X,Y,Z] = ndgrid(-r:r,-r:r,-r:r);
    h = (X.^2/r^2 + Y.^2/r^2 + Z.^2/r^2 <= 1);
    ker = h/sum(h(:));
    imsize = size(mask);
    mask_tmp = convn(mask,ker,'same');
    mask_ero2 = zeros(imsize);
    mask_ero2(mask_tmp > 0.999999) = 1; % no error tolerance
    r = 3; 
    [X,Y,Z] = ndgrid(-r:r,-r:r,-r:r);
    h = (X.^2/r^2 + Y.^2/r^2 + Z.^2/r^2 <= 1);
    ker = h/sum(h(:));
    imsize = size(mask);
    mask_tmp = convn(mask,ker,'same');
    mask_ero3 = zeros(imsize);
    mask_ero3(mask_tmp > 0.999999) = 1; % no error tolerance


    % Mask = mask;
    % Mask_G = Mask;
    % P_B = 30;
    % P = 1 * Mask + P_B * (1-Mask);
    % RDF = 0;
    % save RDF_brain.mat matrix_size voxel_size delta_TE B0_dir CF iMag N_std iFreq Mask Mask_G P RDF
    % % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 500);
    % % nii = make_nii(QSM.*Mask,vox);
    % % save_nii(nii,'TFI/TFI_0_lambda500_full.nii');
    % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1000);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,'TFI/TFI_0_lambda1000_full.nii');
    % % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1500);
    % % nii = make_nii(QSM.*Mask,vox);
    % % save_nii(nii,'TFI/TFI_0_lambda1500_full.nii');
    % % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 2000);
    % % nii = make_nii(QSM.*Mask,vox);
    % % save_nii(nii,'TFI/TFI_0_lambda2000_full.nii');


    % Mask = mask_ero1;
    % Mask_G = Mask;
    % P_B = 30;
    % P = 1 * Mask + P_B * (1-Mask);
    % RDF = 0;
    % save RDF_brain.mat matrix_size voxel_size delta_TE B0_dir CF iMag N_std iFreq Mask Mask_G P RDF
    % % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 500);
    % % nii = make_nii(QSM.*Mask,vox);
    % % save_nii(nii,'TFI/TFI_0_lambda500_ero1.nii');
    % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1000);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,'TFI/TFI_0_lambda1000_ero1.nii');
    % % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1500);
    % % nii = make_nii(QSM.*Mask,vox);
    % % save_nii(nii,'TFI/TFI_0_lambda1500_ero1.nii');
    % % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 2000);
    % % nii = make_nii(QSM.*Mask,vox);
    % % save_nii(nii,'TFI/TFI_0_lambda2000_ero1.nii');


    % Mask = mask_ero2;
    % Mask_G = Mask;
    % P_B = 30;
    % P = 1 * Mask + P_B * (1-Mask);
    % RDF = 0;
    % save RDF_brain.mat matrix_size voxel_size delta_TE B0_dir CF iMag N_std iFreq Mask Mask_G P RDF
    % % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 500);
    % % nii = make_nii(QSM.*Mask,vox);
    % % save_nii(nii,'TFI/TFI_0_lambda500_ero2.nii');
    % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1000);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,'TFI/TFI_0_lambda1000_ero2.nii');
    % % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1500);
    % % nii = make_nii(QSM.*Mask,vox);
    % % save_nii(nii,'TFI/TFI_0_lambda1500_ero2.nii');
    % % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 2000);
    % % nii = make_nii(QSM.*Mask,vox);
    % % save_nii(nii,'TFI/TFI_0_lambda2000_ero2.nii');


    Mask = mask_ero3;
    Mask_G = Mask;
    P_B = 30;
    P = 1 * Mask + P_B * (1-Mask);
    RDF = 0;
    save RDF_brain.mat matrix_size voxel_size delta_TE B0_dir CF iMag N_std iFreq Mask Mask_G P RDF
    % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 500);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,'TFI/TFI_0_lambda500_ero3.nii');
    QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1000);
    nii = make_nii(QSM.*Mask,vox);
    save_nii(nii,'TFI/TFI_0_lambda1000_ero3.nii');
    % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1500);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,'TFI/TFI_0_lambda1500_ero3.nii');
    % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 2000);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,'TFI/TFI_0_lambda2000_ero3.nii');
end




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % complex fitting first then laplacian, resharp, iLSQR
% [iFreq_raw N_std] = Fit_ppm_complex(mag.*exp(1j*ph_corr));
% nii = make_nii(iFreq_raw,vox);
% save_nii(nii,'iFreq_raw.nii');

% delta_TE = TE(2) - TE(1);

% % laplacian unwrapping
% Options.voxelSize = vox;
% iFreq_lap = lapunwrap(iFreq_raw, Options);
% tfs_lap = -iFreq_lap/(2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6);
% nii = make_nii(tfs_lap,vox);
% save_nii(nii,'tfs_lap_ppm_fit.nii');
% % FUDGE unwrapping
% iFreq_fudge = fudge(iFreq_raw);
% tfs_fudge = -iFreq_fudge/(2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6);
% nii = make_nii(tfs_fudge,vox);
% save_nii(nii,'tfs_fudge_ppm_fit.nii');
% % prelude unwrapping
% !prelude -a src/mag1.nii -p iFreq_raw.nii -u iFreq_prelude.nii -m BET_mask.nii -n 12
% !gunzip -f iFreq_prelude.nii.gz
% nii = load_nii(['iFreq_prelude.nii']);
% iFreq_prelude = double(nii.img);
% tfs_prelude = -iFreq_prelude/(2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6);
% nii = make_nii(tfs_prelude,vox);
% save_nii(nii,'tfs_prelude_ppm_fit.nii');
% % bestpath unwrapping
% fid = fopen(['iFreq_raw.dat'],'w');
% fwrite(fid,iFreq_raw,'float');
% fclose(fid);
% disp('--> unwrap aliasing phase using bestpath...');
% mask_unwrp = uint8(abs(mask)*255);
% fid = fopen('mask_unwrp.dat','w');
% fwrite(fid,mask_unwrp,'uchar');
% fclose(fid);
% [pathstr, ~, ~] = fileparts(which('3DSRNCP.m'));
% setenv('pathstr',pathstr);
% setenv('nv',num2str(imsize(1)));
% setenv('np',num2str(imsize(2)));
% setenv('ns',num2str(imsize(3)));
% if isdeployed
%     bash_script = ['~/bin/3DSRNCP iFreq_raw.dat mask_unwrp.dat ' ...
%     'iFreq_bestpath.dat $nv $np $ns reliability${echo_num}.dat'];
% else    
%     bash_script = ['${pathstr}/3DSRNCP iFreq_raw.dat mask_unwrp.dat ' ...
%     'iFreq_bestpath.dat $nv $np $ns reliability${echo_num}.dat'];
% end
% unix(bash_script) ;
% fid = fopen(['iFreq_bestpath.dat'],'r');
% tmp = fread(fid,'float');
% % tmp = tmp - tmp(1);
% iFreq_bestpath = reshape(tmp - round(mean(tmp(mask==1))/(2*pi))*2*pi ,imsize(1:3)).*mask;
% fclose(fid);
% tfs_bestpath = -iFreq_bestpath/(2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6);
% nii = make_nii(tfs_bestpath,vox);
% save_nii(nii,'tfs_bestpath_ppm_fit.nii');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% background field removal

% % (1) lap
% tfs = tfs_lap;
% [lfs_resharp_0, mask_resharp_0] = resharp(tfs,mask,vox,smv_rad,tik_reg,cgs_num);
% % % 3D 2nd order polyfit to remove any residual background
% % lfs_resharp= (lfs_resharp - poly3d(lfs_resharp,mask_resharp)).*mask_resharp;
% % save nifti
% mkdir('RESHARP');
% nii = make_nii(lfs_resharp_0,vox);
% save_nii(nii,['RESHARP/lfs_resharp_0_smvrad' num2str(smv_rad) '_cpx_lap.nii']);
% lfs_resharp_cpx_lap = lfs_resharp_0;

% % (2) FUDGE
% tfs = tfs_fudge;
% [lfs_resharp_0, mask_resharp_0] = resharp(tfs,mask,vox,smv_rad,tik_reg,cgs_num);
% % % 3D 2nd order polyfit to remove any residual background
% % lfs_resharp= (lfs_resharp - poly3d(lfs_resharp,mask_resharp)).*mask_resharp;
% % save nifti
% mkdir('RESHARP');
% nii = make_nii(lfs_resharp_0,vox);
% save_nii(nii,['RESHARP/lfs_resharp_0_smvrad' num2str(smv_rad) '_cpx_fudge.nii']);
% lfs_resharp_cpx_fudge = lfs_resharp_0;


% % (3) prelude
% tfs = tfs_prelude;
% [lfs_resharp_0, mask_resharp_0] = resharp(tfs,mask,vox,smv_rad,tik_reg,cgs_num);
% % % 3D 2nd order polyfit to remove any residual background
% % lfs_resharp= (lfs_resharp - poly3d(lfs_resharp,mask_resharp)).*mask_resharp;
% % save nifti
% mkdir('RESHARP');
% nii = make_nii(lfs_resharp_0,vox);
% save_nii(nii,['RESHARP/lfs_resharp_0_smvrad' num2str(smv_rad) '_cpx_prelude.nii']);
% lfs_resharp_cpx_prelude = lfs_resharp_0;


% % (4) bestpath
% tfs = tfs_bestpath;
% [lfs_resharp_0, mask_resharp_0] = resharp(tfs,mask,vox,smv_rad,tik_reg,cgs_num);
% % % 3D 2nd order polyfit to remove any residual background
% % lfs_resharp= (lfs_resharp - poly3d(lfs_resharp,mask_resharp)).*mask_resharp;
% % save nifti
% mkdir('RESHARP');
% nii = make_nii(lfs_resharp_0,vox);
% save_nii(nii,['RESHARP/lfs_resharp_0_smvrad' num2str(smv_rad) '_cpx_bestpath.nii']);
% lfs_resharp_cpx_bestpath = lfs_resharp_0;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% dipole inversion
% % iLSQR
% % (1) lap
% lfs_resharp_0 = lfs_resharp_cpx_lap;
% chi_iLSQR_0 = QSM_iLSQR(lfs_resharp_0*(2.675e8*dicom_info.MagneticFieldStrength)/1e6,mask_resharp_0,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
% nii = make_nii(chi_iLSQR_0,vox);
% save_nii(nii,['RESHARP/chi_iLSQR_0_niter50_smvrad' num2str(smv_rad) '_cpx_lap.nii']);

% % (2) FUDGE
% lfs_resharp_0 = lfs_resharp_cpx_fudge;
% chi_iLSQR_0 = QSM_iLSQR(lfs_resharp_0*(2.675e8*dicom_info.MagneticFieldStrength)/1e6,mask_resharp_0,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
% nii = make_nii(chi_iLSQR_0,vox);
% save_nii(nii,['RESHARP/chi_iLSQR_0_niter50_smvrad' num2str(smv_rad) '_cpx_fudge.nii']);

% % (3) prelude
% lfs_resharp_0 = lfs_resharp_cpx_prelude;
% chi_iLSQR_0 = QSM_iLSQR(lfs_resharp_0*(2.675e8*dicom_info.MagneticFieldStrength)/1e6,mask_resharp_0,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
% nii = make_nii(chi_iLSQR_0,vox);
% save_nii(nii,['RESHARP/chi_iLSQR_0_niter50_smvrad' num2str(smv_rad) '_cpx_prelude.nii']);

%     % TVDI method
%     sus_resharp = tvdi(lfs_resharp_0,mask_resharp_0,vox,tv_reg,mag(:,:,:,end),z_prjs,inv_num); 
%     nii = make_nii(sus_resharp.*mask_resharp_0,vox);
%     save_nii(nii,'RESHARP/sus_resharp_cpx_prelude.nii');

% % (4) bestpath
% lfs_resharp_0 = lfs_resharp_cpx_bestpath;
% chi_iLSQR_0 = QSM_iLSQR(lfs_resharp_0*(2.675e8*dicom_info.MagneticFieldStrength)/1e6,mask_resharp_0,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
% nii = make_nii(chi_iLSQR_0,vox);
% save_nii(nii,['RESHARP/chi_iLSQR_0_niter50_smvrad' num2str(smv_rad) '_cpx_bestpath.nii']);




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % using only the first echo
% tfs_TE1 = -unph(:,:,:,1)/(2.675e8*dicom_info.MagneticFieldStrength*TE(1)*1e-6);
% [lfs_resharp_0, mask_resharp_0] = -resharp(tfs_TE1,mask,vox,smv_rad,tik_reg,cgs_num);
% % % 3D 2nd order polyfit to remove any residual background
% % lfs_resharp= (lfs_resharp - poly3d(lfs_resharp,mask_resharp)).*mask_resharp;
% % save nifti
% mkdir('RESHARP');
% nii = make_nii(lfs_resharp_0,vox);
% save_nii(nii,['RESHARP/lfs_resharp_0_smvrad' num2str(smv_rad) '_te1.nii']);
% lfs_resharp_te1 = lfs_resharp_0;

% chi_iLSQR_0 = QSM_iLSQR(lfs_resharp_0*(2.675e8*dicom_info.MagneticFieldStrength)/1e6,mask_resharp_0,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
% nii = make_nii(chi_iLSQR_0,vox);
% save_nii(nii,['RESHARP/chi_iLSQR_0_niter50_smvrad' num2str(smv_rad) '_te1.nii']);

% % TVDI method
% sus_resharp = tvdi(lfs_resharp_0,mask_resharp_0,vox,tv_reg,mag(:,:,:,end),z_prjs,inv_num); 
% nii = make_nii(sus_resharp.*mask_resharp,vox);
% save_nii(nii,'RESHARP/sus_resharp_te1.nii');

save('all.mat','-v7.3');
cd(init_dir);

