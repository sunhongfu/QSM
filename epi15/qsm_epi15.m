function qsm_epi15(meas_in, path_out, options)
%QSM_EPI15 Quantitative susceptibility mapping from EPI sequence at 1.5T.
%   QSM_EPI15(MEAS_IN, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   MEAS_IN    - filename or directory of meas file(.out)  : *.out
%   PATH_OUT   - directory to save nifti and/or matrixes   : QSM_EPI_vxxx
%   OPTIONS    - parameter structure including fields below
%    .ph_corr  - N/2 deghosting phase correction method    : 3
%    .ref_coi  - reference coil to use for phase combine   : 3
%    .eig_rad  - radius (mm) of eig decomp kernel          : 4
%    .smv_rad  - radius (mm) of SMV convolution kernel     : 6
%    .tik_reg  - Tikhonov regularization for RESHARP       : 0.001
%    .tv_reg   - Total variation regularization parameter  : 0.0005
%    .bet_thr  - threshold for BET brain mask              : 0.4
%    .tvdi_n   - iteration number of TVDI (nlcg)           : 200
%    .sav_all  - save all the variables for debug          : 0

if ~ exist('meas_in','var') || isempty(meas_in)
    listing = dir([pwd '/*.out']);
    if ~isempty(listing)
        filename = listing(1).name;
        pathstr = pwd;
    else
        error('cannot find meas file');
    end
elseif exist(meas_in,'dir')
    listing = dir([meas_in '/*.out']);
    if ~isempty(listing)
        pathstr = cd(cd(meas_in));
        filename = listing(1).name;
    else
        error('cannot find meas file');
    end
elseif exist(meas_in,'file')
    [pathstr,name,ext] = fileparts(meas_in);
    if isempty(pathstr)
        pathstr = pwd;
    end
    pathstr = cd(cd(pathstr));
    filename = [name ext];
else
    error('cannot find meas file');
end

if ~ exist('path_out','var') || isempty(path_out)
    path_out = pathstr;
end

if ~ exist('options','var') || isempty(options)
    options = [];
end

if ~ isfield(options,'ph_corr')
    options.ph_corr = 3;
    % 1: linear
    % 2: non-linear
    % 3: MTF
end

if ~ isfield(options,'ref_coi')
    options.ref_coi = 3;
end

if ~ isfield(options,'eig_rad')
    options.eig_rad = 4;
end

if ~ isfield(options,'bet_thr')
    options.bet_thr = 0.4;
end

if ~ isfield(options,'smv_rad')
    options.smv_rad = 6;
end

if ~ isfield(options,'tik_reg')
    options.tik_reg = 1e-3;
end

if ~ isfield(options,'tv_reg')
    options.tv_reg = 5e-4;
end

if ~ isfield(options,'tvdi_n')
    options.tvdi_n = 200;
end

if ~ isfield(options,'sav_all')
    options.sav_all = 0;
end

ph_corr = options.ph_corr;
ref_coi = options.ref_coi;
eig_rad = options.eig_rad;
bet_thr = options.bet_thr;
smv_rad = options.smv_rad;
tik_reg = options.tik_reg;
tv_reg  = options.tv_reg;
tvdi_n  = options.tvdi_n;
sav_all = options.sav_all;


% define directories
[~,name] = fileparts(filename);
path_qsm = [path_out, filesep, name '_QSM_EPI15_v200'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);
disp(['Start recon of ' filename]);


% generate raw img
disp('--> reconstruct to complex img ...');
[img,params] = epi15_recon([pathstr,filesep,filename],ph_corr);

% size and resolution
[Nro,Npe,Ns,~] = size(img);
FOV = params.protocol_header.sSliceArray.asSlice{1};
vox = [FOV.dReadoutFOV/Nro, FOV.dPhaseFOV/Npe,  FOV.dThickness];


% % combine coils
cref = 8; % reference coil
radi = 5; % kernel size

img_cmb = sense_se(img,vox,cref,radi);
nii = make_nii(abs(img_cmb),vox);
save_nii(nii,'mag.nii');
nii = make_nii(angle(img_cmb),vox);
save_nii(nii,'ph.nii');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % combine coils
% % 
% img_cmb = zeros(Nro,Npe,Ns);
% matlabpool open
% parfor i = 1:Ns
%     img_cmb(:,:,i) = coilCombinePar(img(:,:,i,:));
% end
% matlabpool close
% nii = make_nii(abs(img_cmb),vox);
% save_nii(nii,'mag.nii');
% nii = make_nii(angle(img_cmb),vox);
% save_nii(nii,'ph.nii');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


! bet mag.nii BET -f 0.4 -m -R;
! gunzip -f BET.nii.gz;
! gunzip -f BET_mask.nii.gz;
nii = load_nii('BET_mask.nii');
mask = double(nii.img);
img_cmb = img_cmb.*mask;

! prelude -a mag.nii -p ph.nii -u unph.nii -m BET_mask.nii -n 8;
! gunzip -f unph.nii.gz;
nii = load_nii('unph.nii');
unph = double(nii.img);


% background field removal
ker_rad = 5; % convolution kernel radius size (mm)
tik_reg = 1e-3; % tikhonov regularization

% % (1) PDF
% theta = -acos(params.protocol_header.sSliceArray.asSlice{1}.sNormal.dTra);
% [lfs,mask_ero] = pdf(tfs,mask,vox,ker_rad,abs(img_cmb),theta);
% nii = make_nii(lfs,vox);
% save_nii(nii,'lfs_pdf.nii');

% (2) RESHARP
[lph,mask_ero] = resharp(unph,mask,vox,ker_rad,tik_reg);
nii = make_nii(mask_ero,vox);
save_nii(nii,'mask_ero.nii');
nii = make_nii(lph,vox);
save_nii(nii,'lph.nii');


% normalize to ppm unit
TE = params.protocol_header.alTE{1}/1e6;
B_0 = params.protocol_header.m_flMagneticFieldStrength;
gamma = 2.675222e8;
lfs = lph/(gamma*TE*B_0)*1e6; % unit ppm
nii = make_nii(lfs,vox);
save_nii(nii,'lfs.nii');


% susceptibility inversion
% account for oblique slicing (head tilted)
theta = -acos(params.protocol_header.sSliceArray.asSlice{1}.sNormal.dTra);

tv_reg = 5e-4; % total variation regularization

sus = tvdi(lfs,mask_ero,vox,tv_reg,abs(img_cmb),theta);
% sus_final = sus.*mask_ero;
% nii = make_nii(sus_final,vox);
nii = make_nii(sus,vox);
save_nii(nii,'sus.nii');


% debugging purpose
% save all the variables in all.mat
save all.mat;


cd(init_dir);