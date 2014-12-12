function qsm_swi15(meas_in, path_out, options)
%QSM_SWI15 Quantitative susceptibility mapping from SWI sequence at 1.5T.
%   QSM_SWI15(MEAS_IN, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   MEAS_IN     - filename or directory of meas file(.out)  : *.out
%   PATH_OUT    - directory to save nifti and/or matrixes   : QSM_*
%   OPTIONS     - parameter structure including fields below
%    .ref_coi   - reference coil to use for phase combine   : 4
%    .eig_rad   - radius (mm) of eig decomp kernel          : 4
%    .ph_unwrap - 'prelude' or 'laplacian' or 'bestpath'    : 'prelude'
%    .smv_rad   - radius (mm) of SMV convolution kernel     : 3
%    .tik_reg   - Tikhonov regularization for RESHARP       : 0.001
%    .tv_reg    - Total variation regularization parameter  : 0.0005
%    .bet_thr   - threshold for BET brain mask              : 0.4
%    .inv_num   - iteration number of TVDI (nlcg)           : 200
%    .save_all  - save all the variables for debug          : 1

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

if ~ isfield(options,'ph_unwrap')
    options.ph_unwrap = 'prelude';
    % % another option is
    % options.ph_unwrap = 'laplacian';
    % % prelude is preferred, unless there's sigularities
    % % in that case, have to use laplacian
    % another option is 'bestpath'
end

if ~ isfield(options,'smv_rad')
    options.smv_rad = 3;
end

if ~ isfield(options,'tik_reg')
    options.tik_reg = 5e-4;
end

if ~ isfield(options,'tv_reg')
    options.tv_reg = 5e-4;
end

if ~ isfield(options,'inv_num')
    options.inv_num = 200;
end

if ~ isfield(options,'save_all')
    options.save_all = 1;
end

if isfield(options,'dicompath')
    dicompath = cd(cd(options.dicompath));
    listing = dir([dicompath, '/*.IMA']);
    dicomfile = [dicompath, '/' listing(1).name];
else
    dicomfile = [];
    setenv('pathstr',pathstr);
    [~,cmout] = unix('find "$pathstr" -name *.IMA | sort');
    if ~ isempty(cmout)
        dicoms = strsplit(cmout,'.IMA');
        dicomfile = [dicoms{1},'.IMA'];
    end
end

ref_coi   = options.ref_coi;
eig_rad   = options.eig_rad;
bet_thr   = options.bet_thr;
ph_unwrap = options.ph_unwrap;
smv_rad   = options.smv_rad;
tik_reg   = options.tik_reg;
tv_reg    = options.tv_reg;
inv_num   = options.inv_num;
save_all  = options.save_all;


% define directories
[~,name] = fileparts(filename);
if strcmpi(ph_unwrap,'prelude')
    path_qsm = [path_out, filesep, 'QSM_SWI15_v500_' name];
elseif strcmpi(ph_unwrap,'laplacian')
    path_qsm = [path_out, filesep, 'QSM_SWI15_v500_lap_' name];
elseif strcmpi(ph_unwrap,'bestpath')
    path_qsm = [path_out, filesep, 'QSM_SWI15_v500_best_' name];
end
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);


% generate raw img
rawfile = {[pathstr,filesep],filename};
[img,params] = swi15_recon(rawfile);


% size and resolution
[nv,np,ns,~] = size(img);
nv = nv;
np = np;
ns = ns;
FOV = params.protocol_header.sSliceArray.asSlice{1};
voxelSize = [FOV.dPhaseFOV/nv, FOV.dReadoutFOV/np,  FOV.dThickness/ns];


% angles!!!
if ~ isempty(dicomfile)
    % read in dicom header, this is accurate information
    info = dicominfo(dicomfile);
    Xz = info.ImageOrientationPatient(3);
    Yz = info.ImageOrientationPatient(6);
    Zz = sqrt(1 - Xz^2 - Yz^2);
    disp('find the dicom');
    dicomfile
    z_prjs = [Xz, Yz, Zz]
else % this would be just an estimation
    sNormal = params.protocol_header.sSliceArray.asSlice{1}.sNormal;
    if ~ isfield(sNormal,'dSag')
        sNormal.dSag = 0;
    end
    if ischar(sNormal.dSag)
        sNormal.dSag = 0;
    end
    if ~ isfield(sNormal,'dCor')
        sNormal.dCor = 0;
    end
    if ischar(sNormal.dCor)
        sNormal.dCor = 0;
    end
    if ~ isfield(sNormal,'dTra')
        sNormal.dTra = 0;
    end
    if ischar(sNormal.dTra)
        sNormal.dTra = 0;
    end
    disp('no dicom found, try to use normal vector');
    z_prjs = [-sNormal.dSag, -sNormal.dCor, sNormal.dTra]
end

% combine RF coils
if size(img,4) > 1
    img_cmb = coils_cmb(img,voxelSize,ref_coi,eig_rad);
else
    img_cmb = img;
end
mkdir('combine');
nii = make_nii(abs(img_cmb),voxelSize);
save_nii(nii,'combine/mag_cmb.nii');
nii = make_nii(angle(img_cmb),voxelSize);
save_nii(nii,'combine/ph_cmb.nii');


% % combine coils (2D SVD)
% img_cmb = zeros([nv,np,ns]);
% matlabpool open
% parfor i = 1:ns
%     img_cmb(:,:,i) = coilCombinePar(img(:,:,i,:));
% end
% matlabpool close
% mkdir('combine');
% nii = make_nii(abs(img_cmb),voxelSize);
% save_nii(nii,'combine/mag_cmb.nii');
% nii = make_nii(angle(img_cmb),voxelSize);
% save_nii(nii,'ph_cmb.nii');


% generate brain mask
setenv('bet_thr',num2str(bet_thr));
unix('bet combine/mag_cmb.nii BET -f ${bet_thr} -m -R');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);



% phase unwrapping, prelude is preferred!
if strcmpi('prelude',ph_unwrap)
    % unwrap combined phase with PRELUDE
    disp('--> unwrap aliasing phase ...');
    bash_script = ['prelude -a combine/mag_cmb.nii -p combine/ph_cmb.nii ' ...
        '-u unph.nii -m BET_mask.nii -n 8'];
    unix(bash_script);
    unix('gunzip -f unph.nii.gz');
    nii = load_nii('unph.nii');
    unph = double(nii.img);


    % % unwrap with Laplacian based method
    % unph = unwrapLaplacian(angle(img_cmb), size(img_cmb), voxelSize);
    % nii = make_nii(unph, voxelSize);
    % save_nii(nii,'unph_lap.nii');


elseif strcmpi('laplacian',ph_unwrap)
    % Ryan Topfer's Laplacian unwrapping
    Options.voxelSize = voxelSize;
    unph = lapunwrap(angle(img_cmb), Options);
    nii = make_nii(unph, voxelSize);
    save_nii(nii,'unph_lap.nii');


elseif strcmpi('bestpath',ph_unwrap)
    % unwrap the phase using best path
    fid = fopen('wrapped_phase.dat','w');
    fwrite(fid,angle(img_cmb),'float');
    fclose(fid);
    % mask_unwrp = uint8(hemo_mask.*255);
    mask_unwrp = uint8(abs(mask)*255);
    fid = fopen('mask_unwrp.dat','w');
    fwrite(fid,mask_unwrp,'uchar');
    fclose(fid);
    unix('cp /home/hongfu/Documents/MATLAB/3DSRNCP 3DSRNCP');
    setenv('nv',num2str(nv));
    setenv('np',num2str(np));
    setenv('ns',num2str(ns));
    bash_script = ['./3DSRNCP wrapped_phase.dat mask_unwrp.dat unwrapped_phase.dat' ...
        '$nv $np $ns reliability.dat'];
    unix(bash_script) ;
    fid = fopen('unwrapped_phase.dat','r');
    unph = fread(fid,'float');
    unph = reshape(unph - unph(1) ,[nv,np,ns]);
    fclose(fid);
    nii = make_nii(unph,voxelSize);
    save_nii(nii,'unph.nii');

    % fid = fopen('reliability.dat','r');
    % reliability = fread(fid,'float');
    % fclose(fid);
    % reliability = reshape(reliability,[nv,np,ns]);
    % reliability = 1./reliability.*mask;
    % reliability_smooth = smooth3(reliability,'box',[7,7,3]);
    % % reliability(reliability <= 0.05) = 0;
    % % reliability(reliability > 0.05) = 1;
    % nii = make_nii(reliability_smooth,voxelSize);
    % save_nii(nii,'reliability_smooth.nii');
else
    error('what unwrapping methods to use? prelude or laplacian or bestpath?')
end


% normalize to ppm unit
TE = params.protocol_header.alTE{1}/1e6;
B_0 = params.protocol_header.m_flMagneticFieldStrength;
gamma = 2.675222e8;
tfs = unph/(gamma*TE*B_0)*1e6; % unit ppm


% background field removal
% (1) RESHARP
[lfs_resharp,mask_resharp] = resharp(tfs,mask,voxelSize,smv_rad,tik_reg);
mkdir('RESHARP');
nii = make_nii(lfs_resharp,voxelSize);
save_nii(nii,'RESHARP/lfs_resharp.nii');

% (2) LBV
lfs_lbv = LBV(tfs,mask,size(unph),voxelSize,0.01,2); % strip 2 layers
mkdir('LBV');
nii = make_nii(lfs_lbv,voxelSize);
save_nii(nii,'LBV/lfs_lbv.nii');
% % Don't use LBV's mask*.bin, not accurate
% % read in eroded mask from LBV
% listing = dir('mask*.bin');
% filename = listing.name;
% fid = fopen(filename);
% mask_lbv = fread(fid,'int');
% mask_lbv = reshape(mask_lbv,size(mask));
% fclose all;
mask_lbv = ones(size(mask));
mask_lbv(lfs_lbv==0) = 0;


% susceptibility inversion
% (1) RESHARP
[sus_resharp,residual_resharp] = tvdi(lfs_resharp,mask_resharp,voxelSize, ...
    tv_reg,abs(img_cmb),z_prjs,inv_num);
nii = make_nii(sus_resharp,voxelSize);
save_nii(nii,'RESHARP/sus_resharp.nii');
nii=make_nii(sus_resharp.*mask_resharp,voxelSize);
save_nii(nii,'RESHARP/sus_resharp_clean.nii');

% (2) LBV
[sus_lbv,residual_lbv] = tvdi(lfs_lbv,mask_lbv,voxelSize,tv_reg, ...
    abs(img_cmb),z_prjs,inv_num);
nii = make_nii(sus_lbv,voxelSize);
save_nii(nii,'LBV/sus_lbv.nii');
nii=make_nii(sus_lbv.*mask_lbv,voxelSize);        
save_nii(nii,'LBV/sus_lbv_clean.nii');  


% save all variables for debugging purpose
if save_all
    clear nii;
    save('all.mat','-v7.3');
end

% save parameters used in the recon
save('parameters.mat','options','-v7.3')


%% clean up
% unix('rm *.nii*');
cd(init_dir);
