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
%    .ref_coi  - reference coil to use for phase combine   : 4
%    .eig_rad  - radius (mm) of eig decomp kernel          : 5
%    .smv_rad  - radius (mm) of SMV convolution kernel     : 6
%    .tik_reg  - Tikhonov regularization for RESHARP       : 0.001
%    .tv_reg   - Total variation regularization parameter  : 0.0005
%    .bet_thr  - threshold for BET brain mask              : 0.4
%    .inv_num  - iteration number of TVDI (nlcg)           : 200
%    .save_all - save all the variables for debug          : 1

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
    options.ref_coi = 4;
end

if ~ isfield(options,'eig_rad')
    options.eig_rad = 4;
end

if ~ isfield(options,'bet_thr')
    options.bet_thr = 0.45;
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

ph_corr  = options.ph_corr;
ref_coi  = options.ref_coi;
eig_rad  = options.eig_rad;
bet_thr  = options.bet_thr;
smv_rad  = options.smv_rad;
tik_reg  = options.tik_reg;
tv_reg   = options.tv_reg;
inv_num  = options.inv_num;
save_all = options.save_all;


% define directories
[~,name] = fileparts(filename);
% path_qsm = [path_out, filesep, strrep(name,' ','_') '_QSM_EPI15_v200'];
path_qsm = [path_out, filesep, 'QSM_EPI15_v100_' name];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);
disp(['Start recon of ' filename]);


% generate raw img
disp('--> reconstruct to complex img ...');
[img,params] = epi15_recon([pathstr,filesep,filename],ph_corr);


% size and resolution
[Nro,Npe,Nsl,~,Nrn] = size(img);
FOV = params.protocol_header.sSliceArray.asSlice{1};
voxelSize = [FOV.dReadoutFOV/Nro, FOV.dPhaseFOV/Npe,  FOV.dThickness];



% angles!!!
if ~ isempty(dicomfile)
    % read in dicom header, this is accurate information
    info = dicominfo(dicomfile);
    Xz = info.ImageOrientationPatient(3);
    Yz = info.ImageOrientationPatient(6);
    Zz = sqrt(1 - Xz^2 - Yz^2);
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
    z_prjs = [-sNormal.dSag, -sNormal.dCor, sNormal.dTra]
end


% save all the variables
if save_all
    img_cmb_all      = zeros([Nro,Npe,Nsl,Nrn]);
    mask_all         = zeros([Nro,Npe,Nsl,Nrn]);
    unph_all         = zeros([Nro,Npe,Nsl,Nrn]);
    lfs_resharp_all  = zeros([Nro,Npe,Nsl,Nrn]);
    mask_resharp_all = zeros([Nro,Npe,Nsl,Nrn]);
    lfs_poly_all     = zeros([Nro,Npe,Nsl,Nrn]);
    sus_resharp_all  = zeros([Nro,Npe,Nsl,Nrn]);
end



% process QSM on individual run volume
img_all = img;
for i = 1:size(img_all,5) % all time series
    img = squeeze(img_all(:,:,:,:,i));

    disp('--> combine multiple channels ...');
    if size(img,4) > 1
        img = coils_cmb(img,voxelSize,ref_coil,eig_rad);
    end

    mkdir('combine');
    nii = make_nii(abs(img),voxelSize);
    save_nii(nii,['combine/mag_cmb' num2str(i,'%03i') '.nii']);
    nii = make_nii(angle(img),voxelSize);
    save_nii(nii,['combine/ph_cmb' num2str(i,'%03i') '.nii']);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % combine coils
    % % 
    % img_cmb = zeros(Nro,Npe,Ns);
    % matlabpool open
    % parfor i = 1:Ns
    %     img_cmb(:,:,i) = coilCombinePar(img(:,:,i,:));
    % end
    % matlabpool close
    % nii = make_nii(abs(img_cmb),voxelSize);
    % save_nii(nii,'mag.nii');
    % nii = make_nii(angle(img_cmb),voxelSize);
    % save_nii(nii,'ph.nii');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    disp('--> extract brain volume and generate mask ...');
    setenv('bet_thr',num2str(bet_thr));
    setenv('time_series',num2str(i,'%03i'));
    [status,cmdout] = unix('rm BET*');
    bash_script = ['bet combine/mag_cmb${time_series}.nii BET${time_series} ' ...
        '-f ${bet_thr} -m -Z'];
    unix(bash_script);
    unix('gunzip -f BET${time_series}.nii.gz');
    unix('gunzip -f BET${time_series}_mask.nii.gz');
    nii = load_nii(['BET' num2str(i,'%03i') '_mask.nii']);
    mask = double(nii.img);

    % unwrap combined phase with PRELUDE
    disp('--> unwrap aliasing phase ...');
    setenv('time_series',num2str(i,'%03i'));
    bash_script = ['prelude -a combine/mag_cmb${time_series}.nii ' ...
        '-p combine/ph_cmb${time_series}.nii -u unph${time_series}.nii ' ...
        '-m BET${time_series}_mask.nii -n 8'];
    unix(bash_script);
    unix('gunzip -f unph${time_series}.nii.gz');
    nii = load_nii(['unph' num2str(i,'%03i') '.nii']);
    unph = double(nii.img);

    % Options.voxelSize = voxelSize;
    % unph = lapunwrap(angle(img), Options);
    % nii = make_nii(unph, voxelSize);
    % save_nii(nii,['unph_lap' num2str(i,'%03i') '.nii']);


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
    save_nii(nii,['RESHARP/lfs_resharp_poly' num2str(i,'%03i') '.nii']);



    disp('--> TV susceptibility inversion ...');
    sus_resharp = tvdi(lfs_resharp,mask_resharp,voxelSize,tv_reg,abs(img), ...
        z_prjs,inv_num);
    nii = make_nii(sus_resharp.*mask_resharp,voxelSize);
    save_nii(nii,['RESHARP/sus_resharp' num2str(i,'%03i') '.nii']);


    % to save all the variables
    if save_all
        img_cmb_all(:,:,:,i)      = img;
        mask_all(:,:,:,i)         = mask;
        unph_all(:,:,:,i)         = unph;
        lfs_resharp_all(:,:,:,i)  = lfs_resharp;
        mask_resharp_all(:,:,:,i) = mask_resharp;
        sus_resharp_all(:,:,:,i)  = sus_resharp;
    end

end

% save all variables for debugging purpose
if save_all
    clear nii;
    save('all.mat','-v7.3');
end

% save parameters used in the recon
save('parameters.mat','options','-v7.3')


% clean up
% unix('rm *.nii*');
cd(init_dir);
