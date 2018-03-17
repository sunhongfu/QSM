function r2sm(path_dicom, path_out, options)
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
%    .tv_reg     - Total variation regularization parameter  : 5e-4
%    .tvdi_n     - iteration number of TVDI (nlcg)           : 200
%    .interp     - interpolate the image to the double size  : 0



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
    options.bkg_rm = {'resharp','lnqsm'};
    % options.bkg_rm = {'pdf','sharp','resharp','esharp','lbv'};
end

if ~ isfield(options,'t_svd')
    options.t_svd = 0.1;
end

if ~ isfield(options,'smv_rad')
    options.smv_rad = 2;
end

if ~ isfield(options,'tik_reg')
    % options.tik_reg = 1e-4;
    options.tik_reg = 0;
end

if ~ isfield(options,'cgs_num')
    options.cgs_num = 500;
end

if ~ isfield(options,'lbv_tol')
    options.lbv_tol = 0.01;
end

if ~ isfield(options,'lbv_peel')
    options.lbv_peel = 1;
end

if ~ isfield(options,'tv_reg')
    options.tv_reg = 5e-4;
end

if ~ isfield(options,'inv_num')
    options.inv_num = 200;
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
list_dicom = list_dicom(~strncmpi('.', {list_dicom.name}, 1));


dicom_info = dicominfo([path_dicom,filesep,list_dicom(1).name]);
dicom_info.EchoTrainLength = 8;

imsize = [dicom_info.Width, dicom_info.Height, ...
            length(list_dicom)/dicom_info.EchoTrainLength/4, ...
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


img = zeros(imsize);
TE = zeros(1, imsize(4));

chopper = double(mod(1:imsize(3),2)) ;
chopper( chopper < 1 ) = -1 ;

Counter = 1;
for zCount = 1 : imsize(3)
    for echoCount = 1 : imsize(4)

		%tmpHeaders{Counter} = dicominfo( imagelist( Counter+2 ).name ) ;
        Counter = Counter + 1 ;
        
        %tmpHeaders{Counter} = dicominfo( imagelist( Counter+2 ).name ) ;
        Counter = Counter + 1 ;
        
        %tmpHeaders{Counter} = dicominfo( imagelist( Counter+2 ).name ) ;
        theReal = ...
            permute(chopper(zCount)*double( dicomread( [path_dicom,filesep,list_dicom(Counter).name] ) ),[2 1]) ;
        dicom_info = dicominfo([path_dicom,filesep,list_dicom(Counter).name]);
	    TE(dicom_info.EchoNumber) = dicom_info.EchoTime*1e-3;
		Counter = Counter + 1 ;
        
        %tmpHeaders{Counter} = dicominfo( imagelist( Counter+2 ).name ) ;
        theImag = ...
            permute(chopper(zCount)*double( dicomread( [path_dicom,filesep,list_dicom(Counter).name] ) ),[2 1]) ;    
        Counter = Counter + 1 ;
        
        img(:,:,zCount,echoCount) = theReal + 1j * theImag ;
    end
end

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

mag = abs(img);
ph = angle(img);
clear img


% define output directories
path_qsm = [path_out '/R2s_SPGR_GE'];
[~,~,~] = mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);

[R2 T2 amp] = r2imgfit(mag, TE);
nii = make_nii(R2,vox);
save_nii(nii,'R2s.nii');
nii = make_nii(T2,vox);
save_nii(nii,'T2s.nii');

