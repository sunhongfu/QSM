% read in uncombined magnitude and phase images
path_mag = '/Users/hongfusun/DATA/7T/1.7.32/1.7.32.30/1.7.32.30.1/1.7.32.30.1.1/dicom_series/25_QSM_gre_9echo_bi_p2_0p75iso';
path_ph = '/Users/hongfusun/DATA/7T/1.7.32/1.7.32.30/1.7.32.30.1/1.7.32.30.1.1/dicom_series/26_QSM_gre_9echo_bi_p2_0p75iso';
path_out = '/Users/hongfusun/DATA/7T/1.7.32/1.7.32.30/1.7.32.30.1/1.7.32.30.1.1/dicom_series';
bet_thr = 0.4;
bet_smooth = 2;
ph_unwrap = 'bestpath';
readout = 'bipolar';
r_mask = 1;
fit_thr = 100;
bkg_rm = 'resharp';
smv_rad = 3;
tik_reg = 1e-6;
cgs_num = 200;
tv_reg = 1e-4;
inv_num = 500;
lbv_tol = 1e-4;
lbv_peel = 3;

% read in DICOMs of both magnitude and raw unfiltered phase images
% read in magnitude DICOMs
path_mag = cd(cd(path_mag));
mag_list = dir(path_mag);
mag_list = mag_list(~strncmpi('.', {mag_list.name}, 1));

% get the sequence parameters
dicom_info = dicominfo([path_mag,filesep,mag_list(1).name]);
EchoTrainLength = 9;
for i = 1:192:1728 % read in TEs
    dicom_info = dicominfo([path_mag,filesep,mag_list(i).name]);
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


mag=zeros(242,280,1728,'single');
for i = 1:length(mag_list)
    mag(:,:,i) = permute(single(dicomread([path_mag,filesep,mag_list(i).name])),[2 1]);
end

% reshape into 9 echoes
mag = reshape(mag,[242,280,192,9]);

% read in magnitude DICOMs
path_ph = cd(cd(path_ph));
ph_list = dir(path_ph);
ph_list = ph_list(~strncmpi('.', {ph_list.name}, 1));


ph=zeros(242,280,1728,'single');
for i = 1:length(ph_list)
    ph(:,:,i) = permute(single(dicomread([path_ph,filesep,ph_list(i).name])),[2 1]);
end

% conver scale to -pi to pi
ph = 2*pi.*ph./4094 - pi;


% reshape into 9 echoes
ph = reshape(ph,[242,280,192,9]);

% save('all.mat','-v7.3');


% size of matrix
imsize = size(mag);


% define output directories
path_qsm = [path_out '/QSM_MEGE_7T_combined'];
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


% generate mask from combined magnitude of the 1th echo
disp('--> extract brain volume and generate mask ...');
setenv('bet_thr',num2str(bet_thr));
setenv('bet_smooth',num2str(bet_smooth));
[status,cmdout] = unix('rm BET*');
unix('bet2 src/mag1.nii BET -f ${bet_thr} -m -w ${bet_smooth}');
% unix('bet2 src/mag1.nii BET -f ${bet_thr} -m');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);


% elseif strcmpi('bipolar',readout)
    ph_corr = zeros(imsize(1:4));
    ph_corr(:,:,:,1:2:end) = geme_cmb(mag(:,:,:,1:2:end,:).*exp(1j*ph(:,:,:,1:2:end,:)),vox,TE(1:2:end),mask);
    ph_corr(:,:,:,2:2:end) = geme_cmb(mag(:,:,:,2:2:end,:).*exp(1j*ph(:,:,:,2:2:end,:)),vox,TE(2:2:end),mask);
% else


% save offset corrected phase niftis
for echo = 1:imsize(4)
    nii = make_nii(ph_corr(:,:,:,echo),vox);
    save_nii(nii,['src/ph_corr' num2str(echo) '.nii']);
end


% % complex fitting without unwrapping of each echo
% iField = mag.*exp(1i*ph_corr);
% %%%%% provide a noise_level here if possible %%%%%%
% if (~exist('noise_level','var'))
%     noise_level = calfieldnoise(iField, mask);
% end
% %%%%% normalize signal intensity by noise to get SNR %%%
% iField = iField/noise_level;
% %%%% Generate the Magnitude image %%%%
% iMag = sqrt(sum(abs(iField).^2,4));
% [iFreq_raw N_std] = Fit_ppm_complex(iField);
% matrix_size = single(imsize(1:3));
% voxel_size = vox;
% delta_TE = TE(2) - TE(1);
% B0_dir = z_prjs';
% CF = dicom_info.ImagingFrequency *1e6;
% tfs = iFreq/(2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6);

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



% fit phase images with echo times
disp('--> magnitude weighted LS fit of phase to TE ...');
<<<<<<< HEAD
[tfs, fit_residual] = echofit(unph,mag_corr,TE,0); 
=======
[tfs, fit_residual] = echofit(unph,mag,TE,0); 
>>>>>>> 2b7f16ba145e12a74d3b6f314aa629fad0d9790b
% [tfs, fit_residual] = echofit(unph,mag,TE,1); 


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
% units: TE s, gamma 2.675e8 rad/(sT), B0 7T
tfs = tfs/(2.675e8*dicom_info.MagneticFieldStrength)*1e6; % unit ppm

nii = make_nii(tfs,vox);
save_nii(nii,'tfs.nii');



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
    lfs_resharp= (lfs_resharp - poly3d(lfs_resharp,mask_resharp)).*mask_resharp;

    % save nifti
    mkdir('RESHARP');
    nii = make_nii(lfs_resharp,vox);
    save_nii(nii,'RESHARP/lfs_resharp.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on RESHARP...');
<<<<<<< HEAD
    sus_resharp = tvdi(lfs_resharp,mask_resharp,vox,tv_reg,mag_corr(:,:,:,end),z_prjs,inv_num); 
   
    % save nifti
    nii = make_nii(sus_resharp.*mask_resharp,vox);
    save_nii(nii,'RESHARP/sus_resharp_1e-5.nii');
=======
    sus_resharp = tvdi(lfs_resharp,mask_resharp,vox,tv_reg,mag(:,:,:,end),z_prjs,inv_num); 
   
    % save nifti
    nii = make_nii(sus_resharp.*mask_resharp,vox);
    save_nii(nii,'RESHARP/sus_resharp.nii');
>>>>>>> 2b7f16ba145e12a74d3b6f314aa629fad0d9790b
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
    lfs_lbv = LBV(tfs,mask.*R,imsize(1:3),vox,lbv_tol,-1,lbv_peel); % strip 2 layers
    mask_lbv = ones(imsize(1:3));
    mask_lbv(lfs_lbv==0) = 0;
    % 3D 2nd order polyfit to remove any residual background
    lfs_lbv = lfs_lbv - poly3d(lfs_lbv,mask_lbv);

    % save nifti
    mkdir('LBV');
    nii = make_nii(lfs_lbv.*mask_lbv,vox);
    save_nii(nii,'LBV/lfs_lbv.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on lbv...');
    sus_lbv = tvdi(lfs_lbv,mask_lbv,vox,tv_reg,mag(:,:,:,end),z_prjs,inv_num);   

    % save nifti
    nii = make_nii(sus_lbv.*mask_lbv,vox);
    save_nii(nii,'LBV/sus_lbv.nii');
end

