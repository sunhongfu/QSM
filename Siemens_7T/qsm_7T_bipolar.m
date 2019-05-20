function qsm_7T_bipolar(path_mag, path_pha, path_out, options)

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
%    .tv_reg     - Total variation regularization parameter  : 5e-4
%    .tvdi_n     - iteration number of TVDI (nlcg)           : 500
%    .interp     - interpolate the image to the double size  : 0


if ~ exist('path_mag','var') || isempty(path_mag)
    error('Please input the directory of magnitude DICOMs')
end

if ~ exist('path_pha','var') || isempty(path_pha)
    error('Please input the directory of unfiltered phase DICOMs')
end

if ~ exist('path_out','var') || isempty(path_out)
    path_out = pwd;
    display('Current directory for output')
end

if ~ exist('options','var') || isempty(options)
    options = [];
end

%% read in DICOMs of both uncombined magnitude and raw unfiltered phase images
path_mag = cd(cd(path_mag));
mag_list = dir([path_mag '/*.dcm']);
mag_list = mag_list(~strncmpi('.', {mag_list.name}, 1));
path_pha = cd(cd(path_pha));
ph_list = dir([path_pha '/*.dcm']);
ph_list = ph_list(~strncmpi('.', {ph_list.name}, 1));

% number of slices (mag and ph should be the same)
nSL = length(ph_list);

% get the sequence parameters
dicom_info = dicominfo([path_pha,filesep,ph_list(end).name]);
NumberOfEchoes = dicom_info.EchoNumber; 

for i = 1:nSL/NumberOfEchoes:nSL % read in TEs
    dicom_info = dicominfo([path_pha,filesep,ph_list(i).name]);
    TE(dicom_info.EchoNumber) = dicom_info.EchoTime*1e-3;
end
vox = [dicom_info.PixelSpacing(1), dicom_info.PixelSpacing(2), dicom_info.SliceThickness];

% angles (z projections of the image x y z coordinates) 
Xz = dicom_info.ImageOrientationPatient(3);
Yz = dicom_info.ImageOrientationPatient(6);
Zxyz = cross(dicom_info.ImageOrientationPatient(1:3),dicom_info.ImageOrientationPatient(4:6));
Zz = Zxyz(3);
z_prjs = [Xz, Yz, Zz];

% read in measurements
mag = zeros(dicom_info.Rows,dicom_info.Columns,nSL,'single');
ph = zeros(dicom_info.Rows,dicom_info.Columns,nSL,'single');
for i = 1:nSL
    mag(:,:,i) = single(dicomread([path_mag,filesep,mag_list(i).name]));
    ph(:,:,i) = single(dicomread([path_pha,filesep,ph_list(i).name]));
end


% crop mosaic into individual images
AcqMatrix = regexp(dicom_info.Private_0051_100b,'(\d)*(\d)','match');
if strcmpi(dicom_info.InPlanePhaseEncodingDirection,'COL') % A/P
% phase encoding along column
    wRow = round(str2num(AcqMatrix{1})/dicom_info.PercentSampling*100);
    wCol = str2num(AcqMatrix{2});
else % L/R
    wCol = round(str2num(AcqMatrix{1})/dicom_info.PercentSampling*100);
    wRow = str2num(AcqMatrix{2});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nCol = double(dicom_info.Columns/wCol);
% nRow = double(dicom_info.Rows/wRow);
% nChan = double(dicom_info.Private_0019_100a);

% mag_all = zeros(wRow,wCol,nChan,nSL,'single');
% ph_all = zeros(wRow,wCol,nChan,nSL,'single');
% for i = 1:nSL
%     for x = 1:wRow
%         for y = 1:wCol
%             for z = 1:nChan
%                 X = floor((z-1)/nCol)*wRow + x;
%                 Y = mod(z-1,nCol)*wCol + y;
%                 mag_all(x,y,z,i) = mag(X,Y,i);
%                 ph_all(x,y,z,i) = ph(X,Y,i);
%             end
%         end
%     end
% end


% % reshape and permute into COLS, ROWS, SLICES, ECHOES, CHANS
% mag_all = reshape(mag_all,[wRow,wCol,nChan,nSL/NumberOfEchoes,NumberOfEchoes]);
% mag_all = permute(mag_all,[2 1 4 5 3]);
% ph_all = reshape(ph_all,[wRow,wCol,nChan,nSL/NumberOfEchoes,NumberOfEchoes]);
% ph_all = permute(ph_all,[2 1 4 5 3]);
% % 0028,0106  Smallest Image Pixel Value: 0
% % 0028,0107  Largest Image Pixel Value: 4094
% % conver scale to -pi to pi
% ph_all = 2*pi.*(ph_all - single(dicom_info.SmallestImagePixelValue))./(single(dicom_info.LargestImagePixelValue - dicom_info.SmallestImagePixelValue)) - pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% works for square 32 channels (faster)
mag = mat2cell(mag,[wRow wRow wRow wRow wRow wRow], [wCol wCol wCol wCol wCol wCol], nSL);
mag_all = cat(4,mag{1,1}, mag{1,2}, mag{1,3}, mag{1,4}, mag{1,5}, mag{1,6}, mag{2,1}, mag{2,2}, mag{2,3}, mag{2,4}, mag{2,5}, mag{2,6}, mag{3,1}, mag{3,2}, mag{3,3}, mag{3,4}, mag{3,5}, mag{3,6}, mag{4,1}, mag{4,2}, mag{4,3}, mag{4,4}, mag{4,5}, mag{4,6}, mag{5,1}, mag{5,2}, mag{5,3}, mag{5,4}, mag{5,5}, mag{5,6}, mag{6,1}, mag{6,2});
clear mag
mag_all = reshape(mag_all, wRow, wCol, nSL/NumberOfEchoes, NumberOfEchoes, 32);
mag_all = permute(mag_all,[2 1 3 4 5]);

ph = mat2cell(ph,[wRow wRow wRow wRow wRow wRow], [wCol wCol wCol wCol wCol wCol], nSL);
ph_all = cat(4,ph{1,1}, ph{1,2}, ph{1,3}, ph{1,4}, ph{1,5}, ph{1,6}, ph{2,1}, ph{2,2}, ph{2,3}, ph{2,4}, ph{2,5}, ph{2,6}, ph{3,1}, ph{3,2}, ph{3,3}, ph{3,4}, ph{3,5}, ph{3,6}, ph{4,1}, ph{4,2}, ph{4,3}, ph{4,4}, ph{4,5}, ph{4,6}, ph{5,1}, ph{5,2}, ph{5,3}, ph{5,4}, ph{5,5}, ph{5,6}, ph{6,1}, ph{6,2});
clear ph
ph_all = reshape(ph_all, wRow, wCol, nSL/NumberOfEchoes, NumberOfEchoes, 32);
ph_all = permute(ph_all,[2 1 3 4 5]);
ph_all = 2*pi.*(ph_all - single(dicom_info.SmallestImagePixelValue))/(single(dicom_info.LargestImagePixelValue - dicom_info.SmallestImagePixelValue)) - pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define output directories
path_qsm = [path_out '/QSM_MEGE_7T'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);
% save the raw data for future use
clear mag ph

imsize = size(ph_all);
save('raw.mat','-v7.3');



% BEGIN THE QSM RECON PIPELINE
% initial quick brain mask
% simple sum-of-square combination
mag1_sos = sqrt(sum(mag_all(:,:,:,1,:).^2,5));
nii = make_nii(mag1_sos,vox);
save_nii(nii,'mag1_sos.nii');

unix('N4BiasFieldCorrection -i mag1_sos.nii -o mag1_sos_n4.nii');

unix('bet2 mag1_sos_n4.nii BET -f 0.2 -m');
% set a lower threshold for postmortem
% unix('bet2 mag1_sos.nii BET -f 0.1 -m');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);

% coil combination % smoothing factor 10?
ph_corr = zeros(imsize(1:4));
mag_corr = zeros(imsize(1:4));

% (1) if unipolar
% [ph_corr,mag_corr] = geme_cmb(mag_all.*exp(1j*ph_all),vox,TE,mask);
% (2) if bipolar
% [ph_corr(:,:,:,1:2:end),mag_corr(:,:,:,1:2:end)] = geme_cmb(mag_all(:,:,:,1:2:end,:).*exp(1j*ph_all(:,:,:,1:2:end,:)),vox,TE(1:2:end),mask,'gaussian');
% !mkdir odd_gaussian; mv offsets* odd_gaussian
% [ph_corr(:,:,:,2:2:end),mag_corr(:,:,:,2:2:end)] = geme_cmb(mag_all(:,:,:,2:2:end,:).*exp(1j*ph_all(:,:,:,2:2:end,:)),vox,TE(2:2:end),mask,'gaussian');
% !mkdir even_gaussian; mv offsets* even_gaussian

[ph_corr(:,:,:,1:2:end),mag_corr(:,:,:,1:2:end)] = geme_cmb(mag_all(:,:,:,1:2:end,:).*exp(1j*ph_all(:,:,:,1:2:end,:)),vox,TE(1:2:end),mask);
!mkdir odd_smooth3; mv offsets* odd_smooth3
[ph_corr(:,:,:,2:2:end),mag_corr(:,:,:,2:2:end)] = geme_cmb(mag_all(:,:,:,2:2:end,:).*exp(1j*ph_all(:,:,:,2:2:end,:)),vox,TE(2:2:end),mask);
!mkdir even_smooth3; mv offsets* even_smooth3

% if need to clear variables for memory
clear mag_all ph_all

% save niftis after coil combination
mkdir('src');
for echo = 1:imsize(4)
    nii = make_nii(mag_corr(:,:,:,echo),vox);
    save_nii(nii,['src/mag_corr' num2str(echo) '.nii']);

    setenv('echo',num2str(echo));
    unix('N4BiasFieldCorrection -i src/mag_corr${echo}.nii -o src/mag_corr${echo}_n4.nii');

    nii = make_nii(ph_corr(:,:,:,echo),vox);
    save_nii(nii,['src/ph_corr' num2str(echo) '.nii']);
end




save('raw.mat','ph_corr','mag_corr','mask','-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) unwrap the phase using best path
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

unph = zeros(imsize(1:4));

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
end

nii = make_nii(unph,vox);
save_nii(nii,'unph_bestpath_before_jump_correction.nii');

% remove all the temp files
! rm *.dat

% 2pi jumps correction
nii = load_nii('unph_diff.nii');
unph_diff = double(nii.img);
unph_diff = unph_diff/2;
for echo = 2:imsize(4)
    meandiff = unph(:,:,:,echo)-unph(:,:,:,1)-double(echo-1)*unph_diff;
    meandiff = meandiff(mask==1);
    meandiff = mean(meandiff(:));
    njump = round(meandiff/(2*pi));
    disp(['    ' num2str(njump) ' 2pi jumps for TE' num2str(echo)]);
    unph(:,:,:,echo) = unph(:,:,:,echo) - njump*2*pi;
    unph(:,:,:,echo) = unph(:,:,:,echo).*mask;
end
nii = make_nii(unph,vox);
save_nii(nii,'unph_bestpath.nii');

unph_bestpath = unph;
save('raw.mat','unph_bestpath','-append');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % (2) unwrap the phase using PRELUDE
% setenv('echo_num',num2str(imsize(4)));
% bash_command = sprintf(['for ph in src/ph_corr[1-$echo_num].nii\n' ...
% 'do\n' ...
% '   base=`basename $ph`;\n' ...
% '   dir=`dirname $ph`;\n' ...
% '   mag=$dir/"mag_corr"${base:7};\n' ...
% '   unph="unph"${base:7};\n' ...
% '   prelude -a $mag -p $ph -u $unph -m BET_mask.nii -n 12&\n' ...
% 'done\n' ...
% 'wait\n' ...
% 'gunzip -f unph*.gz\n']);
% unix(bash_command);

% unph = zeros(imsize(1:4));
% for echo = 1:imsize(4)
%     nii = load_nii(['unph' num2str(echo) '.nii']);
%     unph(:,:,:,echo) = double(nii.img);
% end

% nii = make_nii(unph,vox);
% save_nii(nii,'unph_PRELUDE_before_jump_correction.nii');

% % 2pi jumps correction
% nii = load_nii('unph_diff.nii');
% unph_diff = double(nii.img);
% unph_diff = unph_diff/2;
% for echo = 2:imsize(4)
%     meandiff = unph(:,:,:,echo)-unph(:,:,:,1)-double(echo-1)*unph_diff;
%     meandiff = meandiff(mask==1);
%     meandiff = mean(meandiff(:));
%     njump = round(meandiff/(2*pi));
%     disp(['    ' num2str(njump) ' 2pi jumps for TE' num2str(echo)]);
%     unph(:,:,:,echo) = unph(:,:,:,echo) - njump*2*pi;
%     unph(:,:,:,echo) = unph(:,:,:,echo).*mask;
% end
% nii = make_nii(unph,vox);
% save_nii(nii,'unph_prelude.nii');

% unph_prelude = unph;
% save('raw.mat','unph_prelude','-append');

% mkdir('prelude');
% cd('prelude');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % (3) laplacian unwrapping
% disp('--> unwrap aliasing phase using laplacian...');
% Options.voxelSize = vox;
% for i = 1:imsize(4)
%     unph(:,:,:,i) = lapunwrap(ph_corr(:,:,:,i), Options);
% end
% nii = make_nii(unph, vox);
% save_nii(nii,'unph_lap.nii');

% unph_lap = unph;
% save('raw.mat','unph_lap','-append');

% mkdir('lap');
% cd('lap');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % (4) FUDGE laplacian unwrapping
% disp('--> unwrap aliasing phase using fudge...');
% for i = 1:imsize(4)
%     unph(:,:,:,i) = fudge(ph_corr(:,:,:,i));
% end
% nii = make_nii(unph, vox);
% save_nii(nii,'unph_fudge.nii');

% unph_fudge = unph;
% save('raw.mat','unph_fudge','-append');

% mkdir('fudge');
% cd('fudge');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set parameters
fit_thr = 40;
tik_reg = 1e-6;
cgs_num = 500;
% lsqr_num = 500;
% tv_reg = 2e-4;
% inv_num = 500;


% fit phase images with echo times
disp('--> magnitude weighted LS fit of phase to TE ...');
[tfs_0, fit_residual_0] = echofit(unph,mag_corr,TE,0); 
% normalize to main field
% ph = gamma*dB*TE
% dB/B = ph/(gamma*TE*B0)
% units: TE s, gamma 2.675e8 rad/(sT), B0 7T
tfs_0 = tfs_0/(2.675e8*dicom_info.MagneticFieldStrength)*1e6; % unit ppm
nii = make_nii(tfs_0,vox);
save_nii(nii,'tfs_0.nii');
nii = make_nii(fit_residual_0,vox);
save_nii(nii,'fit_residual_0.nii');

% extra filtering according to fitting residuals
% generate reliability map
fit_residual_0_blur = smooth3(fit_residual_0,'box',round(1./vox)*2+1); 
nii = make_nii(fit_residual_0_blur,vox);
save_nii(nii,'fit_residual_0_blur.nii');
R_0 = ones(size(fit_residual_0_blur));
R_0(fit_residual_0_blur >= fit_thr) = 0;

save('raw.mat','tfs_0','fit_residual_0','fit_residual_0_blur','R_0','fit_thr','cgs_num','tik_reg','-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('raw.mat','tfs_0','mask','R_0','vox','cgs_num','dicom_info','z_prjs','mag_corr','imsize','TE','tik_reg');

for smv_rad = [1 2 3]
    % RE-SHARP (tik_reg: Tikhonov regularization parameter)
    disp('--> RESHARP to remove background field ...');
    % [lfs_resharp_0, mask_resharp_0] = resharp_lsqr(tfs_0,mask.*R_0,vox,smv_rad,lsqr_num);
    [lfs_resharp_0, mask_resharp_0] = resharp(tfs_0,mask.*R_0,vox,smv_rad,tik_reg,cgs_num);
    % save nifti
    mkdir('RESHARP');
    nii = make_nii(lfs_resharp_0,vox);
    % save_nii(nii,['RESHARP/lfs_resharp_0_smvrad' num2str(smv_rad) '_lsqr.nii']);
    save_nii(nii,['RESHARP/lfs_resharp_0_smvrad' num2str(smv_rad) '_cgs_' num2str(tik_reg) '.nii']);

    lfs_resharp(:,:,:,smv_rad) = lfs_resharp_0;
    mask_resharp(:,:,:,smv_rad) = mask_resharp_0;

    % iLSQR
    chi_iLSQR_0 = QSM_iLSQR(lfs_resharp_0*(2.675e8*dicom_info.MagneticFieldStrength)/1e6,mask_resharp_0,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
    nii = make_nii(chi_iLSQR_0,vox);
    save_nii(nii,['RESHARP/chi_iLSQR_smvrad' num2str(smv_rad) '.nii']);

    chi_iLSQR(:,:,:,smv_rad) = chi_iLSQR_0;

    % MEDI
    %%%%% normalize signal intensity by noise to get SNR %%%
    %%%% Generate the Magnitude image %%%%
    iMag = sqrt(sum(mag_corr.^2,4));
    % [iFreq_raw N_std] = Fit_ppm_complex(ph_corr);
    matrix_size = single(imsize(1:3));
    voxel_size = vox;
    delta_TE = TE(2) - TE(1);
    B0_dir = z_prjs';
    CF = dicom_info.ImagingFrequency *1e6;
    iFreq = [];
    N_std = 1;
    RDF = lfs_resharp_0*2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6;
    Mask = mask_resharp_0;
    save RDF.mat RDF iFreq iMag N_std Mask matrix_size...
         voxel_size delta_TE CF B0_dir;
    QSM = MEDI_L1('lambda',1000);
    nii = make_nii(QSM.*Mask,vox);
    save_nii(nii,['RESHARP/MEDI1000_RESHARP_smvrad' num2str(smv_rad) '.nii']);

    chi_MEDI(:,:,:,smv_rad) = QSM.*Mask;
end

save('raw.mat','lfs_resharp','mask_resharp','chi_iLSQR','chi_MEDI','-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for layers = [1 2 3]
%     % LBV
%     disp('--> LBV to remove background field ...');
%     lfs_lbv_0 = LBV(tfs_0,mask.*R_0,imsize(1:3),vox,1e-6,-1,layers); % strip 2 layers
%     mask_lbv_0 = ones(size(mask));
%     mask_lbv_0(lfs_lbv_0==0) = 0;
%     % 3D 2nd order polyfit to remove any residual background
%     lfs_lbv_0= (lfs_lbv_0 - poly3d(lfs_lbv_0,mask_lbv_0)).*mask_lbv_0;

%     % save nifti
%     mkdir('LBV');
%     nii = make_nii(lfs_lbv_0,vox);
%     save_nii(nii,['LBV/lfs_lbv_0_layers' num2str(layers) '.nii']);

%     % iLSQR
%     chi_iLSQR_0 = QSM_iLSQR(lfs_lbv_0*(2.675e8*dicom_info.MagneticFieldStrength)/1e6,mask_resharp_0,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
%     nii = make_nii(chi_iLSQR_0,vox);
%     save_nii(nii,['LBV/chi_iLSQR_0_niter50_layers' num2str(layers) '.nii']);

%     % MEDI
%      %%%%% normalize signal intensity by noise to get SNR %%%
%     %%%% Generate the Magnitude image %%%%
%     iMag = sqrt(sum(mag_corr.^2,4));
%     % [iFreq_raw N_std] = Fit_ppm_complex(ph_corr);
%     matrix_size = single(imsize(1:3));
%     voxel_size = vox;
%     delta_TE = TE(2) - TE(1);
%     B0_dir = z_prjs';
%     CF = dicom_info.ImagingFrequency *1e6;
%     iFreq = [];
%     N_std = 1;
%     RDF = lfs_lbv_0*2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6;
%     Mask = mask_lbv_0;
%     save RDF.mat RDF iFreq iMag N_std Mask matrix_size...
%          voxel_size delta_TE CF B0_dir;
%     QSM = MEDI_L1('lambda',1000);
%     nii = make_nii(QSM.*Mask,vox);
%     save_nii(nii,['LBV/MEDI1000_LBV_layers' num2str(layers) '.nii']);

%     %TVDI
%     sus_resharp = tvdi(lfs_lbv_0,mask_lbv_0,vox,2e-4,iMag,z_prjs,500); 
%     nii = make_nii(sus_resharp.*mask_lbv_0,vox);
%     save_nii(nii,['LBV/TV_2e-4_smvrad' num2str(smv_rad) '.nii']);
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % TFI
% mkdir TFI
% %%%% Generate the Magnitude image %%%%
% iMag = sqrt(sum(mag_corr.^2,4));
% matrix_size = single(imsize(1:3));
% voxel_size = vox;
% delta_TE = TE(2) - TE(1);
% B0_dir = z_prjs';
% CF = dicom_info.ImagingFrequency *1e6;
% N_std = 1;

% % (4) TFI of 3 voxels erosion
% iFreq = tfs_0*2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6;
% % erode the mask (full mask to 3mm erosion)
% % apply R
% mask = mask.*R_0;
% % mask_erosion
% r = 1; 
% [X,Y,Z] = ndgrid(-r:r,-r:r,-r:r);
% h = (X.^2/r^2 + Y.^2/r^2 + Z.^2/r^2 <= 1);
% ker = h/sum(h(:));
% imsize = size(mask);
% mask_tmp = convn(mask,ker,'same');
% mask_ero1 = zeros(imsize);
% mask_ero1(mask_tmp > 0.999999) = 1; % no error tolerance
% r = 2; 
% [X,Y,Z] = ndgrid(-r:r,-r:r,-r:r);
% h = (X.^2/r^2 + Y.^2/r^2 + Z.^2/r^2 <= 1);
% ker = h/sum(h(:));
% imsize = size(mask);
% mask_tmp = convn(mask,ker,'same');
% mask_ero2 = zeros(imsize);
% mask_ero2(mask_tmp > 0.999999) = 1; % no error tolerance
% r = 3; 
% [X,Y,Z] = ndgrid(-r:r,-r:r,-r:r);
% h = (X.^2/r^2 + Y.^2/r^2 + Z.^2/r^2 <= 1);
% ker = h/sum(h(:));
% imsize = size(mask);
% mask_tmp = convn(mask,ker,'same');
% mask_ero3 = zeros(imsize);
% mask_ero3(mask_tmp > 0.999999) = 1; % no error tolerance


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


% Mask = mask_ero3;
% Mask_G = Mask;
% P_B = 30;
% P = 1 * Mask + P_B * (1-Mask);
% RDF = 0;
% save RDF_brain.mat matrix_size voxel_size delta_TE B0_dir CF iMag N_std iFreq Mask Mask_G P RDF
% % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 500);
% % nii = make_nii(QSM.*Mask,vox);
% % save_nii(nii,'TFI/TFI_0_lambda500_ero3.nii');
% QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1000);
% nii = make_nii(QSM.*Mask,vox);
% save_nii(nii,'TFI/TFI_0_lambda1000_ero3.nii');
% % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 1500);
% % nii = make_nii(QSM.*Mask,vox);
% % save_nii(nii,'TFI/TFI_0_lambda1500_ero3.nii');
% % QSM = TFI_L1('filename', 'RDF_brain.mat', 'lambda', 2000);
% % nii = make_nii(QSM.*Mask,vox);
% % save_nii(nii,'TFI/TFI_0_lambda2000_ero3.nii');


% % complex fitting first then laplacian, resharp, iLSQR
% [iFreq_raw N_std] = -Fit_ppm_complex(mag_corr.*exp(1j*ph_corr));
% % laplacian unwrapping
% Options.voxelSize = vox;
% % iFreq = lapunwrap(iFreq_raw, Options);
% iFreq = fudge(iFreq_raw);
% delta_TE = TE(2) - TE(1);
% tfs_0 = iFreq/(2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6);
% smv_rad = 3;
% [lfs_resharp_0, mask_resharp_0] = resharp(tfs_0,mask,vox,smv_rad,tik_reg,cgs_num);
% % % 3D 2nd order polyfit to remove any residual background
% % lfs_resharp= (lfs_resharp - poly3d(lfs_resharp,mask_resharp)).*mask_resharp;
% % save nifti
% mkdir('RESHARP');
% nii = make_nii(lfs_resharp_0,vox);
% save_nii(nii,['RESHARP/lfs_resharp_0_smvrad' num2str(smv_rad) '_cpx_lap.nii']);

% % iLSQR
% chi_iLSQR_0 = QSM_iLSQR(lfs_resharp_0*(2.675e8*dicom_info.MagneticFieldStrength)/1e6,mask_resharp_0,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
% nii = make_nii(chi_iLSQR_0,vox);
% save_nii(nii,['RESHARP/chi_iLSQR_0_niter50_smvrad' num2str(smv_rad) '_cpx_lap.nii']);






% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % perform inversion on each individual echo
% tfs_resharp = zeros(imsize(1:4));
% chi_iLSQR = zeros(imsize(1:4));
% for i = 1:imsize(4)
%     % RE-SHARP (tik_reg: Tikhonov regularization parameter)
%     disp('--> RESHARP to remove background field ...');
%     [tfs_resharp(:,:,:,i), mask_resharp] = resharp(unph(:,:,:,i),mask.*R_0,vox,smv_rad,tik_reg,cgs_num);

%     % iLSQR
%     chi_iLSQR(:,:,:,i) = QSM_iLSQR(tfs_resharp(:,:,:,i),mask_resharp,'H',z_prjs,'voxelsize',vox,'niter',30,'TE',1000*TE(i),'B0',dicom_info.MagneticFieldStrength);
%     nii = make_nii(chi_iLSQR(:,:,:,i),vox);
%     save_nii(nii,['RESHARP/chi_iLSQR_30' 'e' num2str(i) '.nii']);
% end





% % perform inversion on each individual echo
% % use the background from first echo for all echoes
% tfs_resharp = zeros(imsize(1:4));
% chi_iLSQR = zeros(imsize(1:4));
% [tfs_resharp(:,:,:,1), mask_resharp] = resharp(unph(:,:,:,1),mask.*R_0,vox,smv_rad,tik_reg,cgs_num);
% bkg = (unph(:,:,:,1) - tfs_resharp(:,:,:,1)).*mask_resharp;
% for i = 2:imsize(4)
%     tfs_resharp(:,:,:,i) = (unph(:,:,:,i) - bkg/TE(1)*TE(i)).*mask_resharp;
% end
% for i = 1:imsize(4)
%     % RE-SHARP (tik_reg: Tikhonov regularization parameter)
%     % iLSQR
%     chi_iLSQR(:,:,:,i) = QSM_iLSQR(tfs_resharp(:,:,:,i),mask_resharp,'H',z_prjs,'voxelsize',vox,'niter',30,'TE',1000*TE(i),'B0',dicom_info.MagneticFieldStrength);
%     nii = make_nii(chi_iLSQR(:,:,:,i),vox);
%     save_nii(nii,['RESHARP/chi_iLSQR_30_bkg1_' 'e' num2str(i) '.nii']);
% end

