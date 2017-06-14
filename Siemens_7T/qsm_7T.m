% read in uncombined magnitude and phase images
path_mag = '/Users/hongfusun/DATA/7T/1.7-2.32/1.7.32.30/1.7.32.30.1/1.7.32.30.1.1/dicom_series/31_QSM_gre_9echo_bi_p2_0p75iso_QSM_gre_9echo_bi_p2_0p75iso_MAG_UNCOMB';
path_ph = '/Users/hongfusun/DATA/7T/1.7-2.32/1.7.32.30/1.7.32.30.1/1.7.32.30.1.1/dicom_series/32_QSM_gre_9echo_bi_p2_0p75iso_QSM_gre_9echo_bi_p2_0p75iso_PHA_UNCOMB';



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


mag=zeros(1680,1452,1728,'single');
for i = 1:length(mag_list)
    mag(:,:,i) = single(dicomread([path_mag,filesep,mag_list(i).name]));
end



% % mosaic form to 4D nifti
% for i = 1:length(mag_list)
%     mag_mosaic(:,:,i) = single(dicomread([path_mag,filesep,mag_list(i).name]));
% end
% for i = 1:length(ph_list)
%     ph_mosaic(:,:,i) = single(dicomread([path_ph,filesep,ph_list(i).name]));
%     ph_mosaic(:,:,i) = ph_mosaic(:,:,i)/4095*2*pi - pi;
% end

% crop mosaic into individual images
% dicom_info = dicominfo([path_mag,filesep,mag_list(1).name]);
AcqMatrix = regexp(dicom_info.Private_0051_100b,'(\d)*(\d)','match');

if strcmpi(dicom_info.InPlanePhaseEncodingDirection,'COL')
% phase encoding along column
    wRow = round(str2num(AcqMatrix{1})/dicom_info.PercentSampling*100);
    wCol = str2num(AcqMatrix{2});
else
    wCol = round(str2num(AcqMatrix{1})/dicom_info.PercentSampling*100);
    wRow = str2num(AcqMatrix{2});
end

nCol = double(dicom_info.Columns/wCol);
nRow = double(dicom_info.Rows/wRow);
nSL = double(dicom_info.Private_0019_100a);

mag_all = zeros(wRow,wCol,nSL,size(mag,3),'single');
% ph_all = mag_all;
for i = 1:size(mag,3)
    for x = 1:wRow
        for y = 1:wCol
            for z = 1:nSL
                X = floor((z-1)/nCol)*wRow + x;
                Y = mod(z-1,nCol)*wCol + y;
                mag_all(x,y,z,i) = mag(X,Y,i);
                % ph_all(x,y,z,i) = ph_mosaic(X,Y,i);
            end
        end
    end
end

clear mag

% reshape into ROWS, COLS, SLICES, ECHOES, CHANS
mag_all = reshape(mag_all,[wRow,wCol,nSL,192,9]);
mag_all = permute(mag_all,[1 2 4 5 3]);

save('all.mat','-v7.3');


%% phase
clear mag_all

% read in magnitude DICOMs
path_ph = cd(cd(path_ph));
ph_list = dir(path_ph);
ph_list = ph_list(~strncmpi('.', {ph_list.name}, 1));


ph=zeros(1680,1452,1728,'single');
for i = 1:length(ph_list)
    ph(:,:,i) = single(dicomread([path_ph,filesep,ph_list(i).name]));
end

AcqMatrix = regexp(dicom_info.Private_0051_100b,'(\d)*(\d)','match');

if strcmpi(dicom_info.InPlanePhaseEncodingDirection,'COL')
% phase encoding along column
    wRow = round(str2num(AcqMatrix{1})/dicom_info.PercentSampling*100);
    wCol = str2num(AcqMatrix{2});
else
    wCol = round(str2num(AcqMatrix{1})/dicom_info.PercentSampling*100);
    wRow = str2num(AcqMatrix{2});
end

nCol = double(dicom_info.Columns/wCol);
nRow = double(dicom_info.Rows/wRow);
nSL = double(dicom_info.Private_0019_100a);

ph_all = zeros(wRow,wCol,nSL,size(ph,3),'single');
% ph_all = ph_all;
for i = 1:size(ph,3)
    for x = 1:wRow
        for y = 1:wCol
            for z = 1:nSL
                X = floor((z-1)/nCol)*wRow + x;
                Y = mod(z-1,nCol)*wCol + y;
                ph_all(x,y,z,i) = ph(X,Y,i);
                % ph_all(x,y,z,i) = ph_mosaic(X,Y,i);
            end
        end
    end
end

clear ph

% reshape into ROWS, COLS, SLICES, ECHOES, CHANS
ph_all = reshape(ph_all,[wRow,wCol,nSL,192,9]);
ph_all = permute(ph_all,[1 2 4 5 3]);

% 0028,0106  Smallest Image Pixel Value: 0
% 0028,0107  Largest Image Pixel Value: 4094

% conver scale to -pi to pi
ph_all = 2*pi.*ph_all./4094 - pi;

save('all.mat','ph_all','-append');


% size of matrix
imsize = size(mag);


% define output directories
path_qsm = [path_out '/QSM_MEGE_7T'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);




% combine the coils

% combine magnitudes using eig method (DO Walsh, MRM2000)
% if par.nrcvrs > 1
%     disp('--> combine RF rcvrs ...');
%     img_cmb = adaptive_cmb(permute(img,[1 2 3 5 4]),voxelSize,ref_coil,eig_rad);
%     mag_cmb = abs(img_cmb);
%     % at 4.7T, seems the 2rd coil has the best SNR?
% else
%     mag_cmb = abs(img);
% end

% simple sum-of-square combination
mag_all = sum(mag.^2,5);

% save niftis after coil combination
mkdir('src');
for echo = 1:imsize(4)
    nii = make_nii(mag_all(:,:,:,echo),vox);
    save_nii(nii,['src/mag_cmb' num2str(echo) '.nii']);
end


% generate mask from combined magnitude of the 1th echo
disp('--> extract brain volume and generate mask ...');
setenv('bet_thr',num2str(bet_thr));
setenv('bet_smooth',num2str(bet_smooth));
[status,cmdout] = unix('rm BET*');
% unix('bet2 combine/mag_cmb1.nii BET -f ${bet_thr} -m -w ${bet_smooth}');
unix('bet2 src/mag_cmb1.nii BET -f ${bet_thr} -m');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);


% elseif strcmpi('bipolar',readout)
    ph_corr = zeros(imsize(1:4));
    ph_corr(:,:,:,1:2:end) = geme_cmb(mag(:,:,:,1:2:end,:).*exp(1j*ph(:,:,:,1:2:end,:)),vox,TE(1:2:end),mask);
    ph_corr(:,:,:,2:2:end) = geme_cmb(mag(:,:,:,2:2:end,:).*exp(1j*ph(:,:,:,2:2:end,:)),vox,TE(2:2:end),mask);
% else