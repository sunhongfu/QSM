function qsm_spgr_ge_catherine_helix(path_dicom, path_out, corr_echo1)
%QSM_SPGR_GE Quantitative susceptibility mapping from SPGR sequence at GE (3T).
%   QSM_SPGR_GE(PATH_DICOM, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH_DICOM   - directory for input GE dicoms
%   PATH_OUT     - directory to save nifti and/or matrixes   : QSM_SPGR_GE
%   OPTIONS      - parameter structure including fields below
%    .readout    - multi-echo 'unipolar' or 'bipolar'        : 'unipolar'
%    .r_mask     - whether to enable the extra masking       : 1
%    .fit_thr    - extra filtering based on the fit residual : 20
%    .bet_thr    - threshold for BET brain mask              : 0.4
%    .bet_smooth - smoothness of BET brain mask at edges     : 2
%    .ph_unwrap  - 'prelude' or 'bestpath'                   : 'prelude'
%    .bkg_rm     - background field removal method(s)        : 'resharp'
%                  options: 'pdf','sharp','resharp','esharp','lbv'
%                  to try all e.g.: {'pdf','sharp','resharp','esharp','lbv'}
%    .t_svd      - truncation of SVD for SHARP               : 0.1
%    .smv_rad    - radius (mm) of SMV convolution kernel     : 3
%    .tik_reg    - Tikhonov regularization for resharp       : 1e-4
%    .cgs_num    - max interation number for RESHARP         : 200
%    .lbv_peel   - LBV layers to be peeled off               : 2
%    .lbv_tol    - LBV interation error tolerance            : 0.01
%    .tv_reg     - Total variation regularization parameter  : 2e-4
%    .tvdi_n     - iteration number of TVDI (nlcg)           : 500
%    .interp     - interpolate the image to the double size  : 0

if ~ exist('path_dicom','var') || isempty(path_dicom)
    error('Please input the directory of DICOMs')
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
    options.fit_thr = 20;
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
    % options.bkg_rm = {'resharp','lbv','tfi','lnqsm'};
    options.bkg_rm = {'lbv'};
end

if ~ isfield(options,'t_svd')
    options.t_svd = 0.1;
end

if ~ isfield(options,'smv_rad')
    options.smv_rad = 4;
end

if ~ isfield(options,'tik_reg')
    options.tik_reg = 1e-4;
end

if ~ isfield(options,'cgs_num')
    options.cgs_num = 200;
end

if ~ isfield(options,'lbv_tol')
    options.lbv_tol = 0.01;
end

if ~ isfield(options,'lbv_peel')
    options.lbv_peel = 2;
end

if ~ isfield(options,'tv_reg')
    options.tv_reg = 2e-4;
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

% read in MESPGR dicoms (multi-echo gradient-echo)
path_dicom = cd(cd(path_dicom));
list_dicom = dir(path_dicom);

dicom_info = dicominfo([path_dicom,filesep,list_dicom(3).name]);
%dicom_info.EchoTrainLength = 8;

imsize = [dicom_info.Width, dicom_info.Height, (length(list_dicom)-2)/dicom_info.EchoTrainLength/2, ...
			 dicom_info.EchoTrainLength];
vox = [dicom_info.PixelSpacing(1), dicom_info.PixelSpacing(2), dicom_info.SliceThickness];
CF = dicom_info.ImagingFrequency *1e6;

% angles!!!
Xz = dicom_info.ImageOrientationPatient(3);
Yz = dicom_info.ImageOrientationPatient(6);
%Zz = sqrt(1 - Xz^2 - Yz^2);
Zxyz = cross(dicom_info.ImageOrientationPatient(1:3),dicom_info.ImageOrientationPatient(4:6));
Zz = Zxyz(3);
z_prjs = [Xz, Yz, Zz];



Counter = 1;
for zCount = 1 : imsize(3)
    for echoCount = 1 : imsize(4)

        theMag = ...
            permute(double( dicomread( [path_dicom,filesep,list_dicom(Counter+2).name] ) ),[2 1]) ;
        dicom_info = dicominfo([path_dicom,filesep,list_dicom(Counter+2).name]);
	    TE(dicom_info.EchoNumber) = dicom_info.EchoTime*1e-3;
		Counter = Counter + 1 ;
        
        %tmpHeaders{Counter} = dicominfo( imagelist( Counter+2 ).name ) ;
        thePh = ...
            permute(double( dicomread( [path_dicom,filesep,list_dicom(Counter+2).name] ) ),[2 1]) ;    
        Counter = Counter + 1 ;
        
%        img(:,:,zCount,echoCount) = theReal + 1j * theImag ;
    	mag(:,:,zCount,echoCount) = theMag;
    	ph(:,:,zCount,echoCount) = thePh/3141*pi;

	end

end

%%%%%%%%%%%%%%%%%%%%%%%%%
% phase cycling correction
ph(:,:,2:2:end,:)=angle(exp(1j*(ph(:,:,2:2:end,:)+pi)));
%%%%%%%%%%%%%%%%%%%%%%%%%



% define output directories
if strcmpi(corr_echo1,'no')
    path_qsm = [path_out '/QSM_SPGR_GE_no_echo1'];
elseif strcmpi(corr_echo1,'pi')
    path_qsm = [path_out '/QSM_SPGR_GE_pi'];
elseif strcmpi(corr_echo1,'halfpi')
    path_qsm = [path_out '/QSM_SPGR_GE_halfpi'];
end

mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);

% nii = make_nii(mag,vox);
% save_nii(nii,'mag_all.nii');
nii = make_nii(ph,vox);
save_nii(nii,'ph_all_raw.nii');


%%%%%%%%%%%%%%%%%%%%%%%%%
% first echo phase correction
if strcmpi(corr_echo1,'pi')
    ph(:,:,:,1) = angle(exp(1j*(ph(:,:,:,1)-pi)));
    %%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(corr_echo1,'halfpi')
    ph(:,:,:,1) = angle(exp(1j*(ph(:,:,:,1)-pi/2)));
end


% brain extraction
% generate mask from magnitude of the 1th echo
disp('--> extract brain volume and generate mask ...');
setenv('bet_thr',num2str(bet_thr));
setenv('bet_smooth',num2str(bet_smooth));
[status,cmdout] = unix('rm BET*');
nii = make_nii(mag(:,:,:,1),vox);
save_nii(nii,'mag1.nii');
unix('bet2 mag1.nii BET -f ${bet_thr} -m');
% unix('bet2 src/mag5.nii BET -f ${bet_thr} -m -w ${bet_smooth}');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);



% interpolate the images to the double size
if interp
    img = single(img);
    % zero padding the k-space
    k = fftshift(fftshift(fftshift(fft(fft(fft(img,[],1),[],2),[],3),1),2),3);
    k = padarray(k,double(imsize(1:3)/2));
    img = ifft(ifft(ifft(ifftshift(ifftshift(ifftshift(k,1),2),3),[],1),[],2),[],3);
    clear k;
    imsize = size(img);
    vox = vox/2;
end


% save magnitude and raw phase niftis for each echo
mkdir('src')
for echo = 1:imsize(4)
    nii = make_nii(mag(:,:,:,echo),vox);
    save_nii(nii,['src/mag' num2str(echo) '.nii']);
    nii = make_nii(ph(:,:,:,echo),vox);
    save_nii(nii,['src/ph' num2str(echo) '.nii']);
end


if strcmpi(corr_echo1,'no')
    mag = mag(:,:,:,2:end);
    ph = ph(:,:,:,2:end);
    TE = TE(2:end);
    imsize = size(mag);
end


% phase offset correction
% if unipolar
if strcmpi('unipolar',readout)
    ph_corr = geme_cmb(mag.*exp(1j*ph),vox,TE,mask,'gaussian',0);
 % if bipolar
elseif strcmpi('bipolar',readout)
    ph_corr = zeros(imsize);
    ph_corr(:,:,:,1:2:end) = geme_cmb(mag(:,:,:,1:2:end).*exp(1j*ph(:,:,:,1:2:end)),vox,TE(1:2:end),mask,'gaussian',0);
    ph_corr(:,:,:,2:2:end) = geme_cmb(mag(:,:,:,2:2:end).*exp(1j*ph(:,:,:,2:2:end)),vox,TE(2:2:end),mask,'gaussian',0);
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
    '   base=`basename $ph`;\n' ...
    '   dir=`dirname $ph`;\n' ...
    '   mag=$dir/"mag"${base:7};\n' ...
    '   unph="unph"${base:7};\n' ...
    '   prelude -a $mag -p $ph -u $unph -m BET_mask.nii -n 12&\n' ...
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


        if isdeployed
            bash_script = ['~/bin/3DSRNCP wrapped_phase${echo_num}.dat mask_unwrp.dat ' ...
                    'unwrapped_phase${echo_num}.dat $nv $np $ns reliability${echo_num}.dat'];
        else
            bash_script = ['${pathstr}/3DSRNCP wrapped_phase${echo_num}.dat mask_unwrp.dat ' ...
                'unwrapped_phase${echo_num}.dat $nv $np $ns reliability${echo_num}.dat'];
        end
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
if strcmpi('bipolar',readout)
    unph_diff = unph_diff/2;
end 

for echo = 2:imsize(4)
    meandiff = unph(:,:,:,echo)-unph(:,:,:,1)-double(echo-1)*unph_diff;
    meandiff = meandiff(mask==1);
    meandiff = mean(meandiff(:))
    njump = round(meandiff/(2*pi))
    disp(['    ' num2str(njump) ' 2pi jumps for TE' num2str(echo)]);
    unph(:,:,:,echo) = unph(:,:,:,echo) - njump*2*pi;
    unph(:,:,:,echo) = unph(:,:,:,echo).*mask;
end

nii = make_nii(unph,vox);
save_nii(nii,'unph_corr.nii');


% fit phase images with echo times
disp('--> magnitude weighted LS fit of phase to TE ...');
[tfs, fit_residual] = echofit(unph,mag,TE,0); 


% extra filtering according to fitting residuals
if r_mask
    % generate reliability map
    fit_residual_blur = smooth3(fit_residual,'box',round(2./vox)*2+1); 
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
tfs = -tfs/(CF*2*pi)*1e6; % unit ppm

nii = make_nii(tfs,vox);
save_nii(nii,'tfs.nii');


% RE-SHARP (tik_reg: Tikhonov regularization parameter)
if sum(strcmpi('resharp',bkg_rm))

    [lfs_resharp_0, mask_resharp] = resharp(tfs,mask.*R,vox,smv_rad,0,cgs_num);
    % save nifti
    mkdir('RESHARP');
    nii = make_nii(lfs_resharp_0,vox);
    save_nii(nii,['RESHARP/lfs_resharp_tik_', num2str(0), '_num_', num2str(cgs_num), '.nii']);

    % 2D V-SHARP
    voxelsize0 = vox(1:2);
    padsize0 = [100 100]; 
    smvsize = 100;
    for cpt = 1:imsize(3)
      lfs_resharp_v2d(:,:,cpt) = V_SHARP_2d(single(lfs_resharp_0(:,:,cpt)),single(mask_resharp(:,:,cpt)),'smvsize',smvsize,'voxelsize',voxelsize0,'padsize',padsize0);
    end
    mask_resharp_v2d = ones(size(mask));
    mask_resharp_v2d(lfs_resharp_v2d == 0) = 0;

    nii = make_nii(lfs_resharp_v2d,vox);
    save_nii(nii,'RESHARP/lfs_resharp_v2d.nii');

    
    % inversion of susceptibility （iLSQR method)
    chi_iLSQR = QSM_iLSQR(lfs_resharp_v2d*CF*2*pi/1e6,mask_resharp_v2d,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
    nii = make_nii(chi_iLSQR,vox);
    save_nii(nii,['RESHARP/chi_iLSQR_smvrad' num2str(smv_rad) '.nii']);


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
    RDF = lfs_resharp_v2d*CF*2*pi*delta_TE*1e-6;
    Mask = mask_resharp_v2d;
    save RDF.mat RDF iFreq iMag N_std Mask matrix_size...
         voxel_size delta_TE CF B0_dir;
    QSM = MEDI_L1('lambda',1000);
    nii = make_nii(QSM.*Mask,vox);
    save_nii(nii,['RESHARP/MEDI1000_RESHARP_smvrad' num2str(smv_rad) '.nii']);

    % QSM = MEDI_L1('lambda',500);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,['RESHARP/MEDI500_RESHARP_smvrad' num2str(smv_rad) '.nii']);


     iFreq = lfs_resharp_v2d*CF*2*pi*delta_TE*1e-6;
     Mask = mask_resharp_v2d;
     Mask_G = Mask;
     P_B = 30;
     P = 1 * Mask + P_B * (1-Mask);
     RDF = 0;
     save RDF_brain.mat matrix_size voxel_size delta_TE B0_dir CF iMag N_std iFreq Mask Mask_G P RDF
 
     % 1000
     QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1000);
     nii = make_nii(QSM.*Mask,vox);
     save_nii(nii,['RESHARP/TFI_resharp_v2d_lambda1000_smvrad_' num2str(smv_rad) '.nii']);
     
     % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 2000);
     % nii = make_nii(QSM.*Mask,vox);
     % save_nii(nii,'TFI/TFI_lbv_v2d_lambda2000.nii');


    % LN-QSM on local field
    maskR = mask.*R;
    Res_wt = iMag.*maskR;
    Res_wt = Res_wt/sum(Res_wt(:))*sum(maskR(:));

    P = mask_resharp_v2d + 30*(1 - mask_resharp_v2d);
    chi_ero0_500 = tikhonov_qsm(lfs_resharp_v2d, Res_wt.*mask_resharp_v2d, 1, mask_resharp_v2d, mask_resharp_v2d, 0, 1e-4, 0.001, 0, vox, P, z_prjs, 500);
    nii = make_nii(chi_ero0_500.*mask_resharp_v2d,vox);
    save_nii(nii,['RESHARP/chi_ero0_tik_1e-3_tv_1e-4_500_smvrad_' num2str(smv_rad) '.nii']);

end

% LBV
if sum(strcmpi('lbv',bkg_rm))
    disp('--> LBV to remove background field ...');
    mkdir LBV

    % for lbv_peel = 1:2
    for lbv_peel = 1
        lfs_lbv = LBV(tfs,mask.*R,imsize(1:3),vox,lbv_tol,-1,lbv_peel); % strip 1 layer
        mask_lbv = ones(imsize(1:3));
        mask_lbv(lfs_lbv==0) = 0;

        nii = make_nii(lfs_lbv,vox);
        save_nii(nii,'LBV/lfs_lbv.nii');

        voxelsize0 = vox(1:2);
        padsize0 = [100 100]; 
        smvsize = 100;
        for cpt = 1:imsize(3)
          lfs_lbv_v2d(:,:,cpt) = V_SHARP_2d(single(lfs_lbv(:,:,cpt)),single(mask_lbv(:,:,cpt)),'smvsize',smvsize,'voxelsize',voxelsize0,'padsize',padsize0);
        end
        mask_lbv_v2d = ones(size(mask));
        mask_lbv_v2d(lfs_lbv_v2d == 0) = 0;

        nii = make_nii(lfs_lbv_v2d,vox);
        save_nii(nii,'LBV/lfs_lbv_v2d.nii');

      
        % inversion of susceptibility （iLSQR method)
        chi_iLSQR = QSM_iLSQR(lfs_lbv_v2d*(CF*2*pi)/1e6,mask_lbv_v2d,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
        nii = make_nii(chi_iLSQR,vox);
        save_nii(nii,['LBV/chi_iLSQR_peel' num2str(lbv_peel) '.nii']);


    % MEDI
    %%%%% normalize signal intensity by noise to get SNR %%%
    %%%% Generate the Magnitude image %%%%
    iMag = sqrt(sum(mag.^2,4));
    % [iFreq_raw N_std] = Fit_ppm_complex(ph_corr);
    matrix_size = single(imsize(1:3));
    voxel_size = vox;
    delta_TE = TE(2) - TE(1);
    B0_dir = z_prjs';
    iFreq = [];
    N_std = 1;
    RDF = lfs_lbv_v2d*CF*2*pi*delta_TE*1e-6;
    Mask = mask_lbv_v2d;
    save RDF.mat RDF iFreq iMag N_std Mask matrix_size...
         voxel_size delta_TE CF B0_dir;
    % QSM = MEDI_L1('lambda',1000);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,['LBV/MEDI1000_LBV_peel' num2str(lbv_peel) '.nii']);

    QSM = MEDI_L1('lambda',2000);
    nii = make_nii(QSM.*Mask,vox);
    save_nii(nii,['LBV/MEDI2000_LBV_peel' num2str(lbv_peel) '.nii']);


     iFreq = lfs_lbv_v2d*CF*2*pi*delta_TE*1e-6;
     Mask = mask_lbv_v2d;
     Mask_G = Mask;
     P_B = 30;
     P = 1 * Mask + P_B * (1-Mask);
     RDF = 0;
     save RDF_brain.mat matrix_size voxel_size delta_TE B0_dir CF iMag N_std iFreq Mask Mask_G P RDF
 
     % 1000
     % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1000);
     % nii = make_nii(QSM.*Mask,vox);
     % save_nii(nii,['LBV/TFI_lbv_v2d_lambda1000_peel_' num2str(lbv_peel) '.nii']);
     
     % 2000
     QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 2000);
     nii = make_nii(QSM.*Mask,vox);
     save_nii(nii,['LBV/TFI_lbv_v2d_lambda2000_peel_' num2str(lbv_peel) '.nii']);


    % LN-QSM on local field
    maskR = mask.*R;
    Res_wt = iMag.*maskR;
    Res_wt = Res_wt/sum(Res_wt(:))*sum(maskR(:));

    P = mask_lbv_v2d + 30*(1 - mask_lbv_v2d);
    chi_ero0_2000 = tikhonov_qsm(lfs_lbv_v2d, Res_wt.*mask_lbv_v2d, 1, mask_lbv_v2d, mask_lbv_v2d, 0, 1e-4, 0.001, 0, vox, P, z_prjs, 2000);
    nii = make_nii(chi_ero0_2000.*mask_lbv_v2d,vox);
    save_nii(nii,['LBV/chi_ero0_tik_1e-3_tv_1e-4_2000_peel_' num2str(lbv_peel) '.nii']);

    % chi_ero0_2000 = tikhonov_qsm(lfs_lbv_v2d, Res_wt.*mask_lbv_v2d, 1, mask_lbv_v2d, mask_lbv_v2d, 0, 4e-4, 0.001, 0, vox, P, z_prjs, 2000);
    % nii = make_nii(chi_ero0_2000.*mask_lbv_v2d,vox);
    % save_nii(nii,['LBV/chi_ero0_tik_1e-3_tv_4e-4_2000_peel_' num2str(lbv_peel) '.nii']);

    end
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

    %
    iFreq = tfs*CF*2*pi*delta_TE*1e-6;
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


    Mask = mask_ero2;
    Mask_G = Mask;
    P_B = 30;
    P = 1 * Mask + P_B * (1-Mask);
    RDF = 0;
    save RDF_brain.mat matrix_size voxel_size delta_TE B0_dir CF iMag N_std iFreq Mask Mask_G P RDF

    % 1000 looks bad
    QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1000);
    nii = make_nii(QSM.*Mask,vox);
    save_nii(nii,'TFI/TFI_ero2_lambda1000.nii');

    % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 500);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,'TFI/TFI_ero2_lambda500.nii');
end



% LN-QSM
if sum(strcmpi('lnqsm',bkg_rm))
    mkdir LN-QSM

    maskR = mask.*R;
    P = maskR + 30*(1 - maskR);
    Res_wt = iMag.*maskR;
    Res_wt = Res_wt/sum(Res_wt(:))*sum(maskR(:));


    chi_ero0_500 = tikhonov_qsm(tfs, Res_wt.*maskR, 1, maskR, maskR, 0, 1e-4, 0.001, 0, vox, P, z_prjs, 500);
    nii = make_nii(chi_ero0_500,vox);
    save_nii(nii,['LN-QSM/chi_ero0_tik_1e-3_tv_1e-4_500.nii']);


    P = mask_ero1 + 30*(1 - mask_ero1);
    chi_ero1_500 = tikhonov_qsm(tfs, Res_wt.*mask_ero1, 1, mask_ero1, mask_ero1, 0, 1e-4, 0.001, 0, vox, P, z_prjs, 500);
    nii = make_nii(chi_ero1_500,vox);
    save_nii(nii,['LN-QSM/chi_ero1_tik_1e-3_tv_1e-4_500.nii']);


    P = mask_ero2 + 30*(1 - mask_ero2);
    chi_ero2_500 = tikhonov_qsm(tfs, Res_wt.*mask_ero2, 1, mask_ero2, mask_ero2, 0, 1e-4, 0.001, 0, vox, P, z_prjs, 500);
    nii = make_nii(chi_ero2_500,vox);
    save_nii(nii,['LN-QSM/chi_ero2_tik_1e-3_tv_1e-4_500.nii']);


    P = mask_ero3 + 30*(1 - mask_ero3);
    chi_ero3_500 = tikhonov_qsm(tfs, Res_wt.*mask_ero3, 1, mask_ero3, mask_ero3, 0, 1e-4, 0.001, 0, vox, P, z_prjs, 500);
    nii = make_nii(chi_ero3_500,vox);
    save_nii(nii,['LN-QSM/chi_ero3_tik_1e-3_tv_1e-4_500.nii']);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % complex fitting method
% [iFreq_raw N_std] = Fit_ppm_complex(mag.*exp(1j*ph));
% nii = make_nii(iFreq_raw,vox);
% save_nii(nii,'iFreq_raw.nii');
% nii = make_nii(N_std*1e5,vox);
% save_nii(nii,'N_std.nii');

% % % replace laplacian with FUDGE? NO. FUDGE doesn't work
% iFreq_lap = unwrapLaplacian(iFreq_raw, size(iFreq_raw), vox);
% nii = make_nii(iFreq_lap,vox);
% save_nii(nii,'iFreq_lap.nii');

% delta_TE = TE(2) - TE(1);
% tfs = iFreq_lap/(CF*2*pi*delta_TE*1e-6);
% nii = make_nii(tfs,vox);
% save_nii(nii,'tfs_cpx.nii');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % tik-qsm

% % pad zeros
% tfs = padarray(tfs,[0 0 20]);
% mask = padarray(mask,[0 0 20]);
% R = padarray(R,[0 0 20]);

% for r = [1 2 3] 

%     [X,Y,Z] = ndgrid(-r:r,-r:r,-r:r);
%     h = (X.^2/r^2 + Y.^2/r^2 + Z.^2/r^2 <= 1);
%     ker = h/sum(h(:));
%     imsize = size(mask);
%     mask_tmp = convn(mask.*R,ker,'same');
%     mask_ero = zeros(imsize);
%     mask_ero(mask_tmp > 1-1/sum(h(:))) = 1; % no error tolerance

%     % try total field inversion on regular mask, regular prelude
%     Tik_weight = 0.005;
%     TV_weight = 0.003;
%     [chi, res] = tikhonov_qsm(tfs, mask_ero, 1, mask_ero, mask_ero, TV_weight, Tik_weight, vox, z_prjs, 2000);
%     nii = make_nii(chi(:,:,21:end-20).*mask_ero(:,:,21:end-20).*R(:,:,21:end-20),vox);
%     save_nii(nii,['chi_brain_pad20_ero' num2str(r) '_TV_' num2str(TV_weight) '_Tik_' num2str(Tik_weight) '_2000.nii']);

% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% remove files to save space
! rm -r *.dat *.bin *.mat offsets*.nii reliability*.nii results ph_diff.nii unph_diff.nii unph*.nii BET.nii src

% save('all.mat','-v7.3');
cd(init_dir);

