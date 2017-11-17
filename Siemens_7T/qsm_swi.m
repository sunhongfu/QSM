% read in uncombined magnitude and phase images
path_mag = '/home/hongfu/NCIgb5_scratch/hongfu/H315/SWI_H315_MAG_UNCOMB';
path_ph = '/home/hongfu/NCIgb5_scratch/hongfu/H315/SWI_H315_PHA_UNCOMB';
path_out = '/home/hongfu/NCIgb5/hongfu/SWI_H315';

%% read in DICOMs of both uncombined magnitude and raw unfiltered phase images
%
path_mag = cd(cd(path_mag));
mag_list = dir(path_mag);
mag_list = mag_list(~strncmpi('.', {mag_list.name}, 1));
path_ph = cd(cd(path_ph));
ph_list = dir(path_ph);
ph_list = ph_list(~strncmpi('.', {ph_list.name}, 1));

% number of slices (mag and ph should be the same)
nSL = length(ph_list);

% get the sequence parameters
dicom_info = dicominfo([path_ph,filesep,ph_list(1).name]);
TE = dicom_info.EchoTime*1e-3;
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
    ph(:,:,i) = single(dicomread([path_ph,filesep,ph_list(i).name]));
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

nCol = double(dicom_info.Columns/wCol);
nRow = double(dicom_info.Rows/wRow);
nChan = double(dicom_info.Private_0019_100a);

mag_all = zeros(wRow,wCol,nChan,nSL,'single');
ph_all = zeros(wRow,wCol,nChan,nSL,'single');
for i = 1:nSL
    for x = 1:wRow
        for y = 1:wCol
            for z = 1:nChan
                X = floor((z-1)/nCol)*wRow + x;
                Y = mod(z-1,nCol)*wCol + y;
                mag_all(x,y,z,i) = mag(X,Y,i);
                ph_all(x,y,z,i) = ph(X,Y,i);
            end
        end
    end
end

% reshape and permute into COLS, ROWS, SLICES, CHANS
mag_all = reshape(mag_all,[wRow,wCol,nChan,nSL]);
mag_all = permute(mag_all,[2 1 4 3]);
ph_all = reshape(ph_all,[wRow,wCol,nChan,nSL]);
ph_all = permute(ph_all,[2 1 4 3]);
% 0028,0106  Smallest Image Pixel Value: 0
% 0028,0107  Largest Image Pixel Value: 4094
% conver scale to -pi to pi
ph_all = 2*pi.*(ph_all - single(dicom_info.SmallestImagePixelValue))./(single(dicom_info.LargestImagePixelValue - dicom_info.SmallestImagePixelValue)) - pi;

imsize = size(ph_all);

% define output directories
path_qsm = [path_out '/QSM_SWI_7T'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);
% save the raw data for future use
clear mag ph
save('raw.mat','-v7.3');

% adaptive combination
cref = 1;
radi = 5;
poolobj=parpool;
parfor sl = 1:imsize(3)
    img_cmb(:,:,sl) = adaptive_cmb_2d(squeeze(mag_all(:,:,sl,:).*exp(1j*ph_all(:,:,sl,:))),vox,cref,radi);
end
delete(poolobj);


% BEGIN THE QSM RECON PIPELINE
% save niftis after coil combination
mkdir('src');
nii = make_nii(abs(img_cmb),vox);
save_nii(nii,'src/mag_corr.nii');
nii = make_nii(angle(img_cmb),vox);
save_nii(nii,'src/ph_corr.nii');

% brain extraction
unix('bet2 src/mag_corr.nii BET -f 0.1 -m');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);


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

unph = zeros(imsize(1:3));

fid = fopen(['wrapped_phase.dat'],'w');
fwrite(fid,angle(img_cmb),'float');
fclose(fid);

bash_script = ['${pathstr}/3DSRNCP wrapped_phase.dat mask_unwrp.dat ' ...
    'unwrapped_phase.dat $nv $np $ns reliability.dat'];
unix(bash_script) ;

fid = fopen('unwrapped_phase.dat','r');
tmp = fread(fid,'float');
% tmp = tmp - tmp(1);
unph = reshape(tmp - round(mean(tmp(mask==1))/(2*pi))*2*pi ,imsize(1:3)).*mask;
fclose(fid);

nii = make_nii(unph,vox);
save_nii(nii,'unph_bestpath.nii');

% remove all the temp files
! rm *.dat


% Ryan Topfer's Laplacian unwrapping
Options.voxelSize = vox;
unph = lapunwrap(angle(img_cmb), Options).*mask;
nii = make_nii(unph, vox);
save_nii(nii,'unph_lap.nii');


% normalize to main field
% ph = gamma*dB*TE
% dB/B = ph/(gamma*TE*B0)
% units: TE s, gamma 2.675e8 rad/(sT), B0 7T
tfs = unph/(2.675e8*dicom_info.MagneticFieldStrength*TE)*1e6; % unit ppm
nii = make_nii(tfs,vox);
save_nii(nii,'tfs.nii');


% set parameters
smv_rad = 4;
tik_reg = 5e-4;
cgs_num = 50;
tv_reg = 2e-4;
inv_num = 1000;


% RE-SHARP (tik_reg: Tikhonov regularization parameter)
disp('--> RESHARP to remove background field ...');
[lfs_resharp, mask_resharp] = resharp(tfs,mask,vox,smv_rad,tik_reg,cgs_num);
% % 3D 2nd order polyfit to remove any residual background
% lfs_resharp= (lfs_resharp - poly3d(lfs_resharp,mask_resharp)).*mask_resharp;

% save nifti
mkdir('RESHARP');
nii = make_nii(lfs_resharp,vox);
save_nii(nii,'RESHARP/lfs_resharp.nii');

% iLSQR
chi_iLSQR = QSM_iLSQR(lfs_resharp*(2.675e8*dicom_info.MagneticFieldStrength)/1e6,mask_resharp,'H',z_prjs,'voxelsize',vox,'niter',50,'TE',1000,'B0',dicom_info.MagneticFieldStrength);
nii = make_nii(chi_iLSQR,vox);
save_nii(nii,'RESHARP/chi_iLSQR.nii');

% MEDI
%%%%% normalize signal intensity by noise to get SNR %%%
%%%% Generate the Magnitude image %%%%
iMag = abs(img_cmb);
% [iFreq_raw N_std] = Fit_ppm_complex(ph_corr);
matrix_size = single(imsize(1:3));
voxel_size = vox;
delta_TE = 1;
B0_dir = z_prjs';
CF = dicom_info.ImagingFrequency *1e6;
iFreq = [];
N_std = 1;
RDF = lfs_resharp*2.675e8*dicom_info.MagneticFieldStrength*delta_TE*1e-6;
Mask = mask_resharp;
save RDF.mat RDF iFreq iMag N_std Mask matrix_size...
     voxel_size delta_TE CF B0_dir;
% run part of MEDI first
QSM = MEDI_L1('lambda',1000);
nii = make_nii(QSM.*Mask,vox);
save_nii(nii,'RESHARP/MEDI_RESHARP_1000.nii');
QSM = MEDI_L1('lambda',500);
nii = make_nii(QSM.*Mask,vox);
save_nii(nii,'RESHARP/MEDI_RESHARP_500.nii');
