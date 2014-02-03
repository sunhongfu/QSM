function qsm_swi15(meas_in, path_out, options)
%QSM_SWI15 Quantitative susceptibility mapping from SWI sequence at 1.5T.
%   QSM_SWI15(MEAS_IN, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   MEAS_IN    - filename or directory of meas file(.out)  : *.out
%   PATH_OUT   - directory to save nifti and/or matrixes   : QSM_SWI_vxxx
%   OPTIONS    - parameter structure including fields below
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

if ~ isfield(options,'ref_coi')
    options.ref_coi = 4;
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
path_qsm = [path_out, filesep, name '_QSM_SWI15_v200'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);
disp(['Start recon of ' filename]);


% generate raw img
disp('--> reconstruct to complex img ...');
rawfile = {[pathstr,filesep],filename};
[img,params] = swi15_recon(rawfile);


% size and resolution
[Npe,Nro,Ns,~] = size(img);
FOV = params.protocol_header.sSliceArray.asSlice{1};
voxelSize = [FOV.dPhaseFOV/Npe, FOV.dReadoutFOV/Nro,  FOV.dThickness/Ns];


% combine RF coils
disp('--> combine RF rcvrs ...');
img_cmb = sense_se(img,voxelSize,ref_coi,eig_rad);
mkdir('combine');
nii = make_nii(abs(img_cmb),voxelSize);
save_nii(nii,'combine/mag_cmb.nii');
nii = make_nii(angle(img_cmb),voxelSize);
save_nii(nii,'combine/ph_cmb.nii');


% % combine coils (2D SVD)
% img_cmb = zeros([Npe,Nro,Ns]);
% matlabpool open
% parfor i = 1:Ns
%     img_cmb(:,:,i) = coilCombinePar(img(:,:,i,:));
% end
% matlabpool close
% mkdir('combine');
% nii = make_nii(abs(img_cmb),voxelSize);
% save_nii(nii,'combine/mag_cmb.nii');
% nii = make_nii(angle(img_cmb),voxelSize);
% save_nii(nii,'ph_cmb.nii');



% generate brain mask
disp('--> extract brain volume and generate mask ...');
setenv('bet_thr',num2str(bet_thr));
unix('bet combine/mag_cmb.nii BET -f ${bet_thr} -m -R');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);


% unwrap combined phase with PRELUDE
disp('--> unwrap aliasing phase ...');
unix('prelude -a combine/mag_cmb.nii -p combine/ph_cmb.nii -u unph.nii -m BET_mask.nii -n 8');
unix('gunzip -f unph.nii.gz');
nii = load_nii('unph.nii');
unph = double(nii.img);


% background field removal
disp('--> RESHARP to remove background field ...');
mkdir('RESHARP');
[lph_resharp,mask_resharp] = resharp(unph,mask,voxelSize,smv_rad,tik_reg);

% normalize to ppm unit
TE = params.protocol_header.alTE{1}/1e6;
B_0 = params.protocol_header.m_flMagneticFieldStrength;
gamma = 2.675222e8;
lfs_resharp = lph_resharp/(gamma*TE*B_0)*1e6; % unit ppm

nii = make_nii(lfs_resharp,voxelSize);
save_nii(nii,'RESHARP/lfs_resharp.nii');


% susceptibility inversion
disp('--> TV susceptibility inversion ...');
% account for oblique slicing (head tilted)
theta = -acos(params.protocol_header.sSliceArray.asSlice{1}.sNormal.dTra);
sus_resharp = tvdi(lfs_resharp,mask_resharp,voxelSize,tv_reg,abs(img_cmb),theta,tvdi_n);

nii = make_nii(sus_resharp,voxelSize);
save_nii(nii,'RESHARP/sus_resharp.nii');


% save all variables for debugging purpose
if sav_all
    clear nii;
    save('all.mat','-v7.3');
end

% save parameters used in the recon
save('parameters.mat','options','-v7.3')


%% clean up
unix('rm *.nii*');
cd(init_dir);