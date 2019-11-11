function qsm_bravo_uncombined(realDicomsDir, imagDicomsDir, path_out, options, QSM_SPGR_GE_Dir)

num_channel = 12; % number of channels uncombined


if ~ exist('realDicomsDir','var') || isempty(realDicomsDir)
    error('Please input the realDicomsDir')
end

if ~ exist('imagDicomsDir','var') || isempty(imagDicomsDir)
    error('Please input the imagDicomsDir')
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
    options.fit_thr = 5;
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


if ~ isfield(options,'orien')
    options.orien = 'saggital';
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
% interp     = options.interp;



realDicomsDir = cd(cd(realDicomsDir));
real_list = dir(realDicomsDir);
real_list = real_list(~strncmpi('.', {real_list.name}, 1));


% get the sequence parameters
dicom_info = dicominfo([realDicomsDir,filesep,real_list(1).name]);

EchoTrainLength = dicom_info.EchoTrainLength;
for i = 1:length(real_list)./EchoTrainLength/num_channel:length(real_list)./num_channel % read in TEs
    dicom_info = dicominfo([realDicomsDir,filesep,real_list(i).name]);
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


% read in DICOMs
imagDicomsDir = cd(cd(imagDicomsDir));
imag_list = dir(imagDicomsDir);
imag_list = imag_list(~strncmpi('.', {imag_list.name}, 1));

real = zeros(dicom_info.Columns,dicom_info.Rows,length(real_list)./EchoTrainLength/num_channel,EchoTrainLength,num_channel);
imag = zeros(size(real));
for i = 1:length(real_list)
    [NS,NE,NC] = ind2sub([length(real_list)./EchoTrainLength/num_channel,EchoTrainLength,num_channel],i);
    real(:,:,NS,NE,NC) = permute(single(dicomread([realDicomsDir,filesep,real_list(i).name])),[2,1]);
    imag(:,:,NS,NE,NC) = permute(single(dicomread([imagDicomsDir,filesep,imag_list(i).name])),[2,1]);
end

% % correct for pi shift of phase
% imag(:,:,1:2:NS,:,:) = - imag(:,:,1:2:NS,:,:);
% real(:,:,1:2:NS,:,:) = - real(:,:,1:2:NS,:,:);

% size of matrix
imsize = size(real);
if length(imsize)==3
	imsize(4) = 1;
end

% define output directories
path_qsm = [path_out '/QSM_BRAVO_UNCOMBINED'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);

% 
% 
% % interpolate the images to the double size
% if interp
%     img = single(img);
%     % zero padding the k-space
%     k = fftshift(fftshift(fftshift(fft(fft(fft(img,[],1),[],2),[],3),1),2),3);
%     k = padarray(k,double(imsize(1:3)/2));
%     img = ifft(ifft(ifft(ifftshift(ifftshift(ifftshift(k,1),2),3),[],1),[],2),[],3);
%     clear k;
%     imsize = size(img);
%     vox = vox/2;
% end

mag = abs(real+1j*imag);
ph = angle(real+1j*imag);

if strcmpi('unipolar',readout)
    ph = -ph;
end

clear real imag


mag1_sos = squeeze(sqrt(sum(mag(:,:,:,1,:).^2,5)));
nii = make_nii(mag1_sos,vox);
save_nii(nii,'mag1_sos.nii');

iMag = sqrt(sum(squeeze(sqrt(sum(mag.^2,5))).^2,4));
nii = make_nii(iMag,vox);
save_nii(nii,'iMag.nii');

if exist('QSM_SPGR_GE_Dir','var')
    path_qsm = pwd;

    setenv('QSM_SPGR_GE_Dir',QSM_SPGR_GE_Dir);
    setenv('path_qsm',path_qsm);
    
    % register MEGRE to MEBRAVO
    !flirt -in ${QSM_SPGR_GE_Dir}/src/mag1.nii -ref ${path_qsm}/mag1_sos.nii -out ${path_qsm}/qsm2bravo -omat ${path_qsm}/qsm2bravo.mat -bins 256 -cost corratio -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear
    !flirt -in ${QSM_SPGR_GE_Dir}/BET_mask.nii -applyxfm -init ${path_qsm}/qsm2bravo.mat -out ${path_qsm}/BET_mask_flirt -paddingsize 0.0 -interp trilinear -ref ${path_qsm}/mag1_sos.nii
    !gunzip ${path_qsm}/BET_mask_flirt.nii.gz

    % load FLIRT registered BET mask from MEGRE scan
    nii = load_nii('BET_mask_flirt.nii');
    mask = double(nii.img);
    mask = (mask > 0);

else
    % brain extraction
    % generate mask from magnitude of the 1th echo
    disp('--> extract brain volume and generate mask ...');
    setenv('bet_thr',num2str(bet_thr));
    setenv('bet_smooth',num2str(bet_smooth));
    [~,~] = unix('rm BET*');

    % unix('fslswapdim src/mag1.nii -z x -y src/mag1_axial.nii.gz')
    % unix('bet2 src/mag1_axial.nii.gz BET -f 0.5 -m -w 2 -g -0.4');

    if strcmpi('axial',options.orien)
        unix('cp iMag.nii iMag_axial.nii');
        unix('bet2 iMag_axial.nii.gz BET -f 0.3 -m -w 2');
        unix('gunzip -f BET.nii.gz');
        unix('gunzip -f BET_mask.nii.gz');
    else
        unix('fslswapdim iMag.nii -z x -y iMag_axial.nii.gz')
        unix('bet2 iMag_axial.nii.gz BET -f 0.4 -m -w 2 -g -0.4');

        unix('fslswapdim BET_mask.nii.gz y -z -x BET_mask_sag.nii.gz')
        unix('fslswapdim BET.nii.gz y -z -x BET_sag.nii.gz');

        unix('rm BET_mask.nii.gz BET.nii.gz');

        unix('gunzip -f BET_sag.nii.gz');
        unix('gunzip -f BET_mask_sag.nii.gz');

        unix('mv BET_sag.nii BET.nii');
        unix('mv BET_mask_sag.nii BET_mask.nii');
    end

    nii = load_nii('BET_mask.nii');
    mask = double(nii.img);
end


imsize = size(mag);
imsize = imsize(1:4);
ph_corr = zeros(imsize(1:4));
if strcmpi('bipolar',readout)
    ph_corr(:,:,:,1:2:end) = geme_cmb(mag(:,:,:,1:2:end,:).*exp(1j*ph(:,:,:,1:2:end,:)),vox,TE(1:2:end),mask);
    ph_corr(:,:,:,2:2:end) = geme_cmb(mag(:,:,:,2:2:end,:).*exp(1j*ph(:,:,:,2:2:end,:)),vox,TE(2:2:end),mask);
elseif strcmpi('unipolar',readout)
    ph_corr = geme_cmb(mag.*exp(1j*ph),vox,TE,mask);
end


mag_corr = sqrt(sum(mag.^2,5));

clear mag ph


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


elseif strcmpi('laplacian',ph_unwrap)
    % (3) laplacian unwrapping
    disp('--> unwrap aliasing phase using laplacian...');
    Options.voxelSize = vox;
    for i = 1:imsize(4)
        unph(:,:,:,i) = lapunwrap(ph_corr(:,:,:,i), Options).*mask;
    end

    mkdir('lap');
    cd('lap');

    nii = make_nii(unph, vox);
    save_nii(nii,'unph_lap.nii');

    unph_lap = unph;
    % save('raw.mat','unph_lap','-append');



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi('fudge',ph_unwrap)
    % (4) FUDGE laplacian unwrapping
    disp('--> unwrap aliasing phase using fudge...');
    for i = 1:imsize(4)
        unph(:,:,:,i) = fudge(ph_corr(:,:,:,i));
    end
    mkdir('fudge');
    cd('fudge');
    nii = make_nii(unph, vox);
    save_nii(nii,'unph_fudge.nii');

    unph_fudge = unph;
    % save('raw.mat','unph_fudge','-append');


end


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


%%%%%%%%%%%%%%%%%%%%%%%
% insert individual echo recon here
%%%%%%%%%%%%%%%%%%%%%%%



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


% background field removal and dipole inversion
% PDF
if sum(strcmpi('pdf',bkg_rm))
    disp('--> PDF to remove background field ...');
    lfs_pdf = projectionontodipolefields(tfs,mask.*R,vox,mag_corr(:,:,:,end),z_prjs);
    % 3D 2nd order polyfit to remove any residual background
    % lfs_pdf= lfs_pdf - poly3d(lfs_pdf,mask_pdf);

    % save nifti
    [~,~,~] = mkdir('PDF');
    nii = make_nii(lfs_pdf,vox);
    save_nii(nii,'PDF/lfs_pdf.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on PDF...');
    sus_pdf = tvdi(lfs_pdf,mask_pdf,vox,tv_reg,mag_corr(:,:,:,end),z_prjs,inv_num); 

    % save nifti
    nii = make_nii(sus_pdf.*mask_pdf,vox);
    save_nii(nii,['PDF/sus_pdf_tv_', num2str(tv_reg), '_num_', num2str(inv_num), '.nii']);
end

% SHARP (t_svd: truncation threthold for t_svd)
if sum(strcmpi('sharp',bkg_rm))
    disp('--> SHARP to remove background field ...');
    [lfs_sharp, mask_sharp] = sharp(tfs,mask.*R,vox,smv_rad,t_svd);
    % % 3D 2nd order polyfit to remove any residual background
    % lfs_sharp= lfs_sharp - poly3d(lfs_sharp,mask_sharp);

    % save nifti
    [~,~,~] = mkdir('SHARP');
    nii = make_nii(lfs_sharp,vox);
    save_nii(nii,'SHARP/lfs_sharp.nii');
    
    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on SHARP...');
    sus_sharp = tvdi(lfs_sharp,mask_sharp,vox,tv_reg,mag_corr(:,:,:,end),z_prjs,inv_num); 
   
    % save nifti
    nii = make_nii(sus_sharp.*mask_sharp,vox);
    save_nii(nii,['SHARP/sus_sharp_tv_', num2str(tv_reg), '_num_', num2str(inv_num), '.nii']);
end

% RE-SHARP (tik_reg: Tikhonov regularization parameter)
if sum(strcmpi('resharp',bkg_rm))
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
    % sus_resharp = tvdi(lfs_resharp,mask_resharp,vox,tv_reg,mag_corr(:,:,:,end),z_prjs,inv_num); 
   
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
    % QSM = MEDI_L1('lambda',2000);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,['RESHARP/MEDI2000_RESHARP_smvrad' num2str(smv_rad) '.nii']);
    % QSM = MEDI_L1('lambda',1500);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,['RESHARP/MEDI1500_RESHARP_smvrad' num2str(smv_rad) '.nii']);
    % QSM = MEDI_L1('lambda',5000);
    % nii = make_nii(QSM.*Mask,vox);
    % save_nii(nii,['RESHARP/MEDI5000_RESHARP_smvrad' num2str(smv_rad) '.nii']);

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
    % lfs_esharp = lfs_esharp - poly3d(lfs_esharp,mask_esharp);

    % save nifti
    [~,~,~] = mkdir('ESHARP');
    nii = make_nii(lfs_esharp,vox);
    save_nii(nii,'ESHARP/lfs_esharp.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on ESHARP...');
    sus_esharp = tvdi(lfs_esharp,mask_esharp,vox,tv_reg,mag_corr(:,:,:,end),z_prjs,inv_num); 
   
    % save nifti
    nii = make_nii(sus_esharp.*mask_esharp,vox);
    save_nii(nii,['ESHARP/sus_esharp_tv_', num2str(tv_reg), '_num_', num2str(inv_num), '.nii']);
end

% LBV
if sum(strcmpi('lbv',bkg_rm))
   disp('--> LBV to remove background field ...');
   lfs_lbv = LBV(tfs,mask.*R,imsize(1:3),vox,lbv_tol,lbv_peel); % strip 2 layers
   mask_lbv = ones(imsize(1:3));
   mask_lbv(lfs_lbv==0) = 0;
   % 3D 2nd order polyfit to remove any residual background
   lfs_lbv= lfs_lbv - poly3d(lfs_lbv,mask_lbv);

   % save nifti
   [~,~,~] = mkdir('LBV');
   nii = make_nii(lfs_lbv,vox);
   save_nii(nii,'LBV/lfs_lbv.nii');

   % inversion of susceptibility 
   disp('--> TV susceptibility inversion on lbv...');
   sus_lbv = tvdi(lfs_lbv,mask_lbv,vox,tv_reg,mag_corr(:,:,:,end),z_prjs,inv_num);   

   % save nifti
   nii = make_nii(sus_lbv.*mask_lbv,vox);
   save_nii(nii,['LBV/sus_lbv_tv_', num2str(tv_reg), '_num_', num2str(inv_num), '.nii']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tik-qsm
%
%% pad zeros
%tfs_pad = padarray(tfs,[0 0 20]);
%mask_pad = padarray(mask,[0 0 20]);
%R_pad = padarray(R,[0 0 20]);
%
%for r = [1 2 3] 
%
%%    [X,Y,Z] = ndgrid(-r:r,-r:r,-r:r);
%    h = (X.^2/r^2 + Y.^2/r^2 + Z.^2/r^2 <= 1);
%    ker = h/sum(h(:));
%    imsize = size(mask_pad);
%    mask_tmp = convn(mask_pad.*R_pad,ker,'same');
%    mask_ero = zeros(imsize);
%    mask_ero(mask_tmp > 1-1/sum(h(:))) = 1; % no error tolerance
%
%    % try total field inversion on regular mask, regular prelude
%    Tik_weight = 0.008;
%%    TV_weight = 0.003;
%    chi = tikhonov_qsm(tfs_pad, mask_ero, 1, mask_ero, mask_ero, TV_weight, Tik_weight, vox, z_prjs, 2000);
%    nii = make_nii(chi(:,:,21:end-20).*mask_ero(:,:,21:end-20).*R_pad(:,:,21:end-20),vox);
%    save_nii(nii,['chi_brain_pad20_ero' num2str(r) '_TV_' num2str(TV_weight) '_Tik_' num2str(Tik_weight) '_2000.nii']);
%
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% R2* fitting
[R2, T2, amp] = r2imgfit2(double(mag_corr),TE,repmat(mask,[1 1 1 imsize(4)]));
nii = make_nii(R2,vox);
save_nii(nii,'R2.nii');
nii = make_nii(T2,vox);
save_nii(nii,'T2.nii');
nii = make_nii(amp,vox);
save_nii(nii,'amp.nii');

save('all.mat','-v7.3');
cd(init_dir);
