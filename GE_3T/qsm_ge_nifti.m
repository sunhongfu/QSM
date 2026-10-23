function qsm_ge_nifti(path_converted, path_mask, path_out, options)
%QSM_GE_NIFTI Conventional multi-echo GRE QSM from NIfTIs converted by dicom_to_nifti.py.
%   QSM_GE_NIFTI(PATH_CONVERTED, PATH_MASK, PATH_OUT, OPTIONS)
%
%   Pipeline: phase offset correction (geme_cmb) -> best-path unwrapping per echo
%   -> 2pi jump correction -> magnitude-weighted echo fit -> RESHARP
%   -> dipole inversion with TV (tvdi) and, if STI Suite is on the path, iLSQR.
%
%   PATH_CONVERTED - folder written by dicom_to_nifti.py (from the iQSM / iQSM_Plus repos):
%                    dcm_converted_phase_4d.nii.gz, dcm_converted_magnitude_4d.nii.gz, params.json
%   PATH_MASK      - brain mask NIfTI on the same grid (e.g. bet2 on the first-echo magnitude)
%   PATH_OUT       - output directory                          : PATH_CONVERTED/QSM_GE_NIFTI
%   OPTIONS        - parameter structure including fields below
%    .readout    - multi-echo 'unipolar' or 'bipolar'        : 'unipolar'
%    .phase_sign - +1 or -1; GE phase is opposite to Siemens : -1
%                  (same as --reverse-phase-sign 1 in iQSM / iQSM+)
%    .r_mask     - whether to enable the extra masking       : 1
%    .fit_thr    - extra filtering based on the fit residual : 20
%    .smv_rad    - radius (mm) of SMV convolution kernel     : 3
%    .tik_reg    - Tikhonov regularization for resharp       : 1e-4
%    .cgs_num    - max interation number for RESHARP         : 200
%    .tv_reg     - Total variation regularization parameter  : 5e-4
%    .inv_num    - iteration number of TVDI (nlcg)           : 500
%    .ilsqr      - also run iLSQR (needs STI Suite on path)  : 1
%
%   Outputs (NIfTI, ppm, same geometry as the input):
%    tfs_ppm, lfs_resharp_ppm, mask_resharp, QSM_RESHARP_TV, QSM_RESHARP_iLSQR
%
%   Requires MATLAB Image Processing Toolbox (niftiread/niftiwrite) and the platform's
%   3DSRNCP binary in phase_unwrapping/ (see phase_unwrapping/3DSRNCP.m).

if ~ exist('path_out','var') || isempty(path_out)
    path_out = fullfile(path_converted, 'QSM_GE_NIFTI');
end

if ~ exist('options','var') || isempty(options)
    options = [];
end

if ~ isfield(options,'readout')
    options.readout = 'unipolar';
end

if ~ isfield(options,'phase_sign')
    options.phase_sign = -1;
end

if ~ isfield(options,'r_mask')
    options.r_mask = 1;
end

if ~ isfield(options,'fit_thr')
    options.fit_thr = 20;
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

if ~ isfield(options,'tv_reg')
    options.tv_reg = 5e-4;
end

if ~ isfield(options,'inv_num')
    options.inv_num = 500;
end

if ~ isfield(options,'ilsqr')
    options.ilsqr = 1;
end

readout    = options.readout;
phase_sign = options.phase_sign;
r_mask     = options.r_mask;
fit_thr    = options.fit_thr;
smv_rad    = options.smv_rad;
tik_reg    = options.tik_reg;
cgs_num    = options.cgs_num;
tv_reg     = options.tv_reg;
inv_num    = options.inv_num;

path_converted = cd(cd(path_converted));
path_mask = abspath(path_mask);
[~,~,~] = mkdir(path_out);
path_out = cd(cd(path_out));

% read converted NIfTIs and acquisition parameters
p = jsondecode(fileread(fullfile(path_converted, 'params.json')));
ph   = double(niftiread(fullfile(path_converted, 'dcm_converted_phase_4d.nii.gz')));
mag  = double(niftiread(fullfile(path_converted, 'dcm_converted_magnitude_4d.nii.gz')));
mask = double(niftiread(path_mask) > 0);
TE     = p.te_ms(:)' * 1e-3;
vox    = p.voxel_size_mm(:)';
B0     = p.b0_T;
z_prjs = p.b0_direction(:)';
imsize = size(ph);

% output header: 3D, single precision, same geometry as the input
info = niftiinfo(fullfile(path_converted, 'dcm_converted_phase_4d.nii.gz'));
info.ImageSize = imsize(1:3);
info.PixelDimensions = info.PixelDimensions(1:3);
info.raw.dim(1) = 3;
info.raw.dim(5) = 1;
info.Datatype = 'single';
info.BitsPerPixel = 32;
info.Description = 'qsm_ge_nifti';

init_dir = pwd;
cd(path_out);
wr = @(img, name) niftiwrite(single(img), name, info, 'Compressed', true);


% phase offset correction (writes unph_diff.nii in the current folder)
disp('--> phase offset correction ...');
if strcmpi('unipolar',readout)
    ph_corr = geme_cmb(mag.*exp(1j*ph),vox,TE,mask);
elseif strcmpi('bipolar',readout)
    ph_corr = zeros(imsize);
    ph_corr(:,:,:,1:2:end) = geme_cmb(mag(:,:,:,1:2:end).*exp(1j*ph(:,:,:,1:2:end)),vox,TE(1:2:end),mask);
    ph_corr(:,:,:,2:2:end) = geme_cmb(mag(:,:,:,2:2:end).*exp(1j*ph(:,:,:,2:2:end)),vox,TE(2:2:end),mask);
else
    error('is the sequence unipolar or bipolar readout?')
end


% unwrap the phase of each echo using best path
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
    unix(bash_script);

    fid = fopen(['unwrapped_phase' num2str(echo_num) '.dat'],'r');
    tmp = fread(fid,'float');
    unph(:,:,:,echo_num) = reshape(tmp - round(mean(tmp(mask==1))/(2*pi))*2*pi ,imsize(1:3)).*mask;
    fclose(fid);
end


% check and correct for 2pi jump between echoes
disp('--> correct for potential 2pi jumps between TEs ...')
unph_diff = double(niftiread('unph_diff.nii'));
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
[tfs, fit_residual] = echofit(unph,mag,TE,0);


% extra filtering according to fitting residuals
if r_mask
    fit_residual_blur = smooth3(fit_residual,'box',round(1./vox)*2+1);
    R = ones(size(fit_residual_blur));
    R(fit_residual_blur >= fit_thr) = 0;
else
    R = 1;
end


% normalize to main field
% ph = gamma*dB*TE
% dB/B = ph/(gamma*TE*B0)
% units: TE s, gamma 2.675e8 rad/(sT), B0 in T
tfs = phase_sign*tfs/(2.675e8*B0)*1e6; % unit ppm
wr(tfs.*mask, 'tfs_ppm');


% RE-SHARP (tik_reg: Tikhonov regularization parameter)
disp('--> RESHARP to remove background field ...');
[lfs_resharp, mask_resharp] = resharp(tfs,mask.*R,vox,smv_rad,tik_reg,cgs_num);
wr(lfs_resharp, 'lfs_resharp_ppm');
wr(mask_resharp, 'mask_resharp');


% TV inversion of susceptibility
disp('--> TV susceptibility inversion on RESHARP...');
sus_resharp = tvdi(lfs_resharp,mask_resharp,vox,tv_reg,mag(:,:,:,end),z_prjs,inv_num);
wr(sus_resharp.*mask_resharp, 'QSM_RESHARP_TV');


% iLSQR (STI Suite, Li et al. NeuroImage 2015)
if options.ilsqr
    if exist('QSM_iLSQR','file')
        disp('--> iLSQR susceptibility inversion on RESHARP...');
        chi_iLSQR = QSM_iLSQR(lfs_resharp*(2.675e8*B0)/1e6,mask_resharp,'H',z_prjs, ...
            'voxelsize',vox,'niter',50,'TE',1000,'B0',B0);
        wr(chi_iLSQR.*mask_resharp, 'QSM_RESHARP_iLSQR');
    else
        warning('QSM_iLSQR not found -- add STI Suite to the MATLAB path to get the iLSQR map.');
    end
end


[~,~] = unix('rm -f *.dat');
cd(init_dir);
disp(['DONE. Results in ' path_out]);
end


function f = abspath(f)
% absolute path, so the file can still be found after cd(path_out)
if ~ (startsWith(f, filesep) || ~isempty(regexp(f, '^[A-Za-z]:', 'once')))
    f = fullfile(pwd, f);
end
end
