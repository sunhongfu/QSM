nii =  load_nii('recon/unaliasedImage_mag.nii');
mag_all = single(permute(nii.img,[2 1 3 4 5]));

nii =  load_nii('recon/unaliasedImage_pha.nii');
ph_all = single(permute(nii.img,[2 1 3 4 5]));
ph_all(:,:,:,1:2:end,:) = angle(-exp(1j*ph_all(:,:,:,1:2:end,:)));

path_dicom = '/Users/uqhsun8/Downloads/HONGFU_PRODUCTION_1/Subject1/A_A/HbaDevel - 17100/BRAVO_0';

if ~ exist('path_dicom','var') || isempty(path_dicom)
    error('Please input the directory of DICOMs')
end

if ~ exist('path_out','var') || isempty(path_out)
    path_out = pwd;
    disp('Current directory for output')
end

if ~ exist('options','var') || isempty(options)
    options = [];
end

if ~ isfield(options,'readout')
    options.readout = 'bipolar';
end

if ~ isfield(options,'r_mask')
    options.r_mask = 1;
end

if ~ isfield(options,'fit_thr')
    options.fit_thr = 20;
end

if ~ isfield(options,'bet_thr')
    options.bet_thr = 0.3;
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


path_dicom = cd(cd(path_dicom));
real_list = dir(path_dicom);
real_list = real_list(~strncmpi('.', {real_list.name}, 1));

% get the sequence parameters
dicom_info = dicominfo([path_dicom,filesep,real_list(1).name]);

%EchoTrainLength = dicom_info.EchoTrainLength;
EchoTrainLength =4;
for i = 1:EchoTrainLength % read in TEs
    dicom_info = dicominfo([path_dicom,filesep,real_list(i).name]);
    TE(dicom_info.EchoNumber) = dicom_info.EchoTime*1e-3;
end
vox = [dicom_info.PixelSpacing(1), dicom_info.PixelSpacing(2), dicom_info.SliceThickness];


% angles!!! (z projections)
Xz = dicom_info.ImageOrientationPatient(3);
Yz = dicom_info.ImageOrientationPatient(6);
%Zz = sqrt(1 - Xz^2 - Yz^2);
Zxyz = cross(dicom_info.ImageOrientationPatient(1:3),dicom_info.ImageOrientationPatient(4:6));
Zz = Zxyz(3);
z_prjs = [Xz, Yz, Zz];


% define output directories
path_qsm = [path_out '/QSM_BRAVO'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);


iMag = squeeze(sqrt(sum(mag_all(:,:,:,:,1).^2,3)));
nii = make_nii(iMag,vox);
save_nii(nii,'iMag.nii');




% brain extraction
% generate mask from magnitude of the 1th echo
disp('--> extract brain volume and generate mask ...');
setenv('bet_thr',num2str(bet_thr));
setenv('bet_smooth',num2str(bet_smooth));
[~,~] = unix('rm BET*');
unix('bet2 iMag.nii BET -f 0.4 -m -w 2');
%unix('bet2 src/mag1.nii BET -f 0.4 -m -w 2');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);

mag_all = permute(mag_all,[1 2 4 5 3]);
ph_all = permute(ph_all,[1 2 4 5 3]);
imsize = size(mag_all);
imsize = imsize(1:4)
ph_corr = zeros(imsize(1:4));
ph_corr(:,:,:,1:2:end) = geme_cmb(mag_all(:,:,:,1:2:end,:).*exp(1j*ph_all(:,:,:,1:2:end,:)),vox,TE(1:2:end),mask);
ph_corr(:,:,:,2:2:end) = geme_cmb(mag_all(:,:,:,2:2:end,:).*exp(1j*ph_all(:,:,:,2:2:end,:)),vox,TE(2:2:end),mask);

mag_corr = sqrt(sum(mag_all.^2,5));

% save offset corrected phase niftis
mkdir('src')
for echo = 1:imsize(4)
    nii = make_nii(ph_corr(:,:,:,echo),vox);
    save_nii(nii,['src/ph_corr' num2str(echo) '.nii']);
end
for echo = 1:imsize(4)
    nii = make_nii(mag_corr(:,:,:,echo),vox);
    save_nii(nii,['src/mag_corr' num2str(echo) '.nii']);
end


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

imsize = imsize(1:4)
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



if imsize(4) > 1
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
        meandiff = mean(meandiff(:));
        njump = round(meandiff/(2*pi));
        disp(['    ' num2str(njump) ' 2pi jumps for TE' num2str(echo)]);
        unph(:,:,:,echo) = unph(:,:,:,echo) - njump*2*pi;
        unph(:,:,:,echo) = unph(:,:,:,echo).*mask;
    end

    % for echo = 2:imsize(4)
    %     meandiff = unph(:,:,:,echo)-unph(:,:,:,echo-1)-unph_diff;
    %     meandiff = meandiff(mask==1);
    %     meandiff = mean(meandiff(:))
    %     njump = round(meandiff/(2*pi))
    %     disp(['    ' num2str(njump) ' 2pi jumps for TE' num2str(echo)]);
    %     unph(:,:,:,echo) = unph(:,:,:,echo) - njump*2*pi;
    %     unph(:,:,:,echo) = unph(:,:,:,echo).*mask;
    % end

    nii = make_nii(unph,vox);
    save_nii(nii,'unph_corrected.nii');
end



% fit phase images with echo times
disp('--> magnitude weighted LS fit of phase to TE ...');
if imsize(4) > 1
	[tfs, fit_residual] = echofit(unph,mag_corr,TE,0); 
else
	tfs = unph/TE;
end

if imsize(4) > 1
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
else
    R = 1;
end




% normalize to main field
% ph = gamma*dB*TE
% dB/B = ph/(gamma*TE*B0)
% units: TE s, gamma 2.675e8 rad/(sT), B0 3T
if (strcmp(readout,'unipolar'))
    tfs = tfs/(2.675e8*dicom_info.MagneticFieldStrength)*1e6; % unit ppm
else 
    tfs = -tfs/(2.675e8*dicom_info.MagneticFieldStrength)*1e6; % unit ppm
end

nii = make_nii(tfs,vox);
save_nii(nii,'tfs.nii');



disp('--> RESHARP to remove background field ...');
    [lfs_resharp, mask_resharp] = resharp(tfs,mask.*R,vox,smv_rad,tik_reg,cgs_num);
    % % 3D 2nd order polyfit to remove any residual background
    % lfs_resharp= lfs_resharp - poly3d(lfs_resharp,mask_resharp);

    % save nifti
    [~,~,~] = mkdir('RESHARP');
    nii = make_nii(lfs_resharp,vox);
    save_nii(nii,['RESHARP/lfs_resharp_tik_', num2str(tik_reg), '_num_', num2str(cgs_num), '.nii']);

    % % inversion of susceptibility 
    % disp('--> TV susceptibility inversion on RESHARP...');
    % sus_resharp = tvdi(lfs_resharp,mask_resharp,vox,tv_reg,mag(:,:,:,end),z_prjs,inv_num); 
   
    % % save nifti
    % nii = make_nii(sus_resharp.*mask_resharp,vox);
    % save_nii(nii,['RESHARP/sus_resharp_tik_', num2str(tik_reg), '_tv_', num2str(tv_reg), '_num_', num2str(inv_num), '.nii']);
    
    
    % iLSQR
    chi_iLSQR = QSM_iLSQR(lfs_resharp*(2.675e8*dicom_info.MagneticFieldStrength)/1e6,mask_resharp,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
    nii = make_nii(chi_iLSQR,vox);
    save_nii(nii,['RESHARP/chi_iLSQR_smvrad' num2str(smv_rad) '.nii']);
    
    % MEDI
    %%%%% normalize signal intensity by noise to get SNR %%%
    %%%% Generate the Magnitude image %%%%
    if imsize(4) > 1
        iMag = sqrt(sum(mag_corr.^2,4));
	else
		iMag = mag_corr;
    end
    
    % [iFreq_raw N_std] = Fit_ppm_complex(ph_corr);
    matrix_size = single(imsize(1:3));
    voxel_size = vox;
    if imsize(4) > 1
        delta_TE = TE(2) - TE(1);
    else
        delta_TE = 0.001;
    end
    B0_dir = z_prjs';
    CF = dicom_info.ImagingFrequency *1e6;
    iFreq = [];
    N_std = 1;
    RDF = lfs_resharp*2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6;
    Mask = mask_resharp;
    save RDF.mat RDF iFreq iMag N_std Mask matrix_size...
         voxel_size delta_TE CF B0_dir;
    QSM = MEDI_L1('lambda',1000);
    nii = make_nii(QSM.*Mask,vox);
    save_nii(nii,['RESHARP/MEDI1000_RESHARP_smvrad' num2str(smv_rad) '.nii']);
    QSM = MEDI_L1('lambda',2000);
    nii = make_nii(QSM.*Mask,vox);
    save_nii(nii,['RESHARP/MEDI2000_RESHARP_smvrad' num2str(smv_rad) '.nii']);
    QSM = MEDI_L1('lambda',1500);
    nii = make_nii(QSM.*Mask,vox);
    save_nii(nii,['RESHARP/MEDI1500_RESHARP_smvrad' num2str(smv_rad) '.nii']);
    QSM = MEDI_L1('lambda',5000);
    nii = make_nii(QSM.*Mask,vox);
    save_nii(nii,['RESHARP/MEDI5000_RESHARP_smvrad' num2str(smv_rad) '.nii']);

    save('all.mat','-v7.3');
    