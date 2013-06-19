function qsm_r2s47(path_in, path_out, params)
%QSM_R2S47 Quantitative susceptibility mapping from R2* sequence at 4.7T.
%   QSM_R2S47(PATH_IN, PATH_OUT, PARAMS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH_IN    - directory of .fid from gemsme3d sequence  : pwd/gemsme3d*.fid
%   PATH_OUT   - directory to save nifti and/or matrixes   : gemsme3d*.fid/QSM
%   PARAMS     - parameter structure including fields below (!in small case!)
%    .save_mat - whether to save matrixes (1) or not (0)   : 0
%    .bkgrm    - background field removal method(s)        : {'resharp','pdf'}
%    .ker_rad  - radius (mm) of RESHARP convolution kernel : 6
%    .tik_reg  - Tikhonov regularization for RESHARP       : 0.001
%    .tsvd     - truncation of SVD for SHARP
%    .tv_reg   - Total variation regularization parameter  : 0.0005


%% default settings and prompts
if ~ exist('path_in','var') || isempty(path_in)
    if exist('fid','file')
        path_in = pwd;
    elseif exist('gemsme3d_R2s_01.fid','dir')
        path_in = [pwd, '/gemsme3d_R2s_01.fid'];
    else
        error('define fid directory.');
    end
end

if ~ exist('path_out','var') || isempty(path_out)
    path_out = [path_in '/QSM'];
end

if ~ exist('params','var') || isempty(params)
    params = [];
end

if ~ isfield(params,'bkgrm')
    params.bkgrm = {'sharp','resharp','pdf'};
end

if ~ isfield(params,'ker_rad')
    params.ker_rad = 6;
end

if ~ isfield(params,'tik_reg')
    params.tik_reg = 1e-3;
end

if ~ isfield(params,'tsvd')
    params.tsvd = 0.05;
end

if ~ isfield(params,'tv_reg')
    params.tv_reg = 5e-4;
end

if ~ isfield(params,'save_mat')
    params.save_mat = 0;
end

bkgrm    = params.bkgrm;
ker_rad  = params.ker_rad;
tik_reg  = params.tik_reg;
tsvd     = params.tsvd;
tv_reg   = params.tv_reg;
save_mat = params.save_mat;


%% define directories
mkdir(path_out);

if save_mat
    path_mat = [path_out '/matrix'];
    mkdir(path_mat);
end

path_nft = [path_out '/nifti'];
mkdir(path_nft);

TEMP = [path_out '/temp'];
mkdir(TEMP);

disp('... switch to the temp folder ...');
cd(TEMP);


%% reconstruct complex image from fid file
disp('--> (1/9) reconstruct fid to complex img ...');
[img,par] = reconfid(path_in);

% save matrix
if save_mat
    mkdir([path_mat '/rawdata']);
    save([path_mat '/rawdata/img.mat'],'img','-v7.3');
    save([path_mat '/rawdata/par.mat'],'par','-v7.3');
end

% keep only the first 4 echoes (last 15.2ms)
img = img(:,:,:,1:4,:);
par.ne = 4;
[np,nv,nv2,ne,~] = size(img);
voxelSize = par.res; % resolution in mm/pixel


%% combine magnitude channels using marc's 'arrayrec.m'
disp('--> (2/9) combine rcvrs for magnitude ...');
mag_cmb = zeros(np,nv,nv2,ne);
for echo = 1:ne
    mag_cmb(:,:,:,echo) = arrayrec(squeeze(img(:,:,:,echo,:)),1/2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DO: remove coil sense profile (Hammond paper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save matrix
if save_mat
    mkdir([path_mat '/combine']);
    save([path_mat '/combine/mag_cmb.mat'],'mag_cmb','-v7.3');
end

% save nifti
mkdir([path_nft '/combine']);
for echo =  1:ne
    nii = make_nii(mag_cmb(:,:,:,echo),voxelSize);
    save_nii(nii,[path_nft '/combine/mag_te' num2str(echo) '.nii']);
end


%% generate mask from SOS combined magnitude of the last echo
disp('--> (3/9) extract brain volume and generate mask ...');
setenv('path_nft', path_nft);
unix('bet $path_nft/combine/mag_te4.nii BET -f 0.3 -m -R');
%unix('bet $path_nft/combine/mag_te1.nii BET -f 0.5 -m -R');
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii('BET_mask.nii');
mask = double(nii.img);

% save matrix
if save_mat
    mkdir([path_mat '/mask']);
    save([path_mat '/mask/mask.mat'],'mask','-v7.3');
end

% save nifti
mkdir([path_nft '/mask']);
copyfile('BET_mask.nii',[path_nft '/mask/mask.nii']);


%% combine phase channels
disp('--> (4/9) combine rcvrs for phase ...');
% (1) filter smooth the complex offsets
ph_cmb = sense(img,par);
% (2) unwrap the offsets and filter with medfilt3
%ph_cmb = sense_med(img,par);

% save matrix
if save_mat
    save([path_mat '/combine/ph_cmb.mat'],'ph_cmb','-v7.3');
end

% save nifti
for echo =  1:ne
    nii = make_nii(ph_cmb(:,:,:,echo),voxelSize);
    save_nii(nii,[path_nft '/combine/ph_te' num2str(echo) '.nii']);
end

clear img;


%% unwrap phase from each echo
disp('--> (5/9) unwrap aliasing phase for 4 TEs ...');

bash_command = sprintf(['for ph in $path_nft/combine/ph*\n' ...
'do\n' ...
'	base=`basename $ph`;\n' ...
'	dir=`dirname $ph`;\n' ...
'	mag=$dir/"mag"${base:2};\n' ...
'	unph="unph"${base:2};\n' ...
'	prelude -a $mag -p $ph -u $unph -m BET_mask.nii -n 8&\n' ...
'done\n' ...
'wait\n' ...
'gunzip -f unph*.gz\n']);

unix(bash_command);

unph_cmb = zeros(np,nv,nv2,ne);
for echo = 1:ne
    nii = load_nii(['unph_te' num2str(echo) '.nii']);
    unph_cmb(:,:,:,echo) = double(nii.img);
end


% check and correct for 2pi jump between echoes
disp('--> (6/9) correct for potential 2pi jumps between TEs ...')

nii = load_nii('unph_diff.nii');
unph_diff = double(nii.img);

for echo = 2:ne
    meandiff = unph_cmb(:,:,:,echo)-unph_cmb(:,:,:,1)-(echo-1)*unph_diff;
    meandiff = meandiff(mask==1);
    meandiff = mean(meandiff(:));
    njump = round(meandiff/(2*pi));
    disp(['    ' num2str(njump) ' 2pi jumps for TE' num2str(echo)]);
    unph_cmb(:,:,:,echo) = unph_cmb(:,:,:,echo) - njump*2*pi;
end

% save matrix
if save_mat
    mkdir([path_mat '/unwrap']);
    save([path_mat '/unwrap/unph_cmb.mat'],'unph_cmb','-v7.3'); 
end

% save nifti
mkdir([path_nft '/unwrap']);
for echo = 1:ne
    nii = make_nii(unph_cmb(:,:,:,echo),voxelSize);
    save_nii(nii,[path_nft '/unwrap/unph_te' num2str(echo) '.nii']);
end

clear ph_cmb


%% fit phase images with echo times
disp('--> (7/9) magnitude weighted LS fit of phase to TE ...');
[tfs, fit_residual] = echofit(unph_cmb,mag_cmb,par); 

% normalize to main field
% ph = gamma*dB*TE
% dB/B = ph/(gamma*TE*B0)
% units: TE s, gamma 2.675e8 rad/(sT), B0 4.7T
tfs = -tfs/(2.675e8*4.7)*1e6; % unit ppm

% generate reliability map
R = ones(size(fit_residual));
R(fit_residual >= 5) = 0;

% save matrix
if save_mat
    mkdir([path_mat '/fit']);
    save([path_mat '/fit/tfs.mat'],'tfs','-v7.3');
    save([path_mat '/fit/fit_residual.mat'],'fit_residual','-v7.3');
    save([path_mat '/fit/R.mat'],'R','-v7.3');
end

% save nifti
mkdir([path_nft '/fit']);
nii = make_nii(tfs,voxelSize);
save_nii(nii,[path_nft '/fit/tfs.nii']);
nii = make_nii(fit_residual,voxelSize);
save_nii(nii,[path_nft '/fit/fit_residual.nii']);
nii = make_nii(R,voxelSize);
save_nii(nii,[path_nft '/fit/R.nii']);

clear unph_cmb


%% SHARP (tsvd: truncation threthold for TSVD)
if sum(strcmpi('sharp',bkgrm))
    disp('--> (8/9) SHARP to remove background field ...');
    [lfs, mask_ero] = sharp(tfs,mask.*R,voxelSize,ker_rad,tsvd);
    mask_final = mask_ero;

    % save matrix
    if save_mat
        mkdir([path_mat '/rmbkg']);
        save([path_mat '/rmbkg/lfs_sharp.mat'],'lfs','-v7.3');
        save([path_mat '/mask/mask_sharp_final.mat'],'mask_final','-v7.3');
    end

    % save nifti
    mkdir([path_nft '/rmbkg/']);
    nii = make_nii(lfs,voxelSize);
    save_nii(nii,[path_nft '/rmbkg/lfs_sharp_xy.nii']);
    nii = make_nii(permute(lfs,[1 3 2]),voxelSize);
    save_nii(nii,[path_nft '/rmbkg/lfs_sharp_xz.nii']);
    nii = make_nii(permute(lfs,[2 3 1]),voxelSize);
    save_nii(nii,[path_nft '/rmbkg/lfs_sharp_yz.nii']);
    nii = make_nii(mask_final,voxelSize);
    save_nii(nii,[path_nft '/mask/mask_sharp_final.nii']);

    

    % inversion of susceptibility 
    disp('--> (9/9) TV susceptibility inversion on SHARP...');
    sus = tvdi(lfs, mask_final, voxelSize, tv_reg, mag_cmb(:,:,:,4)); 
   
    % save matrix
    if save_mat
        mkdir([path_mat '/inversion']);
        save([path_mat '/inversion/sus_sharp.mat'],'sus','-v7.3');
    end

    % save nifti
    mkdir([path_nft '/inversion']);
    nii = make_nii(sus,voxelSize);
    save_nii(nii,[path_nft '/inversion/sus_sharp_xy.nii']);
    nii = make_nii(permute(sus,[1 3 2]),voxelSize);
    save_nii(nii,[path_nft '/inversion/sus_sharp_xz.nii']);
    nii = make_nii(permute(sus,[2 3 1]),voxelSize);
    save_nii(nii,[path_nft '/inversion/sus_sharp_yz.nii']);
end


%% RE-SHARP (tik_reg: Tikhonov regularization parameter)
if sum(strcmpi('resharp',bkgrm))
    disp('--> (8/9) RE-SHARP to remove background field ...');
    [lfs, mask_ero] = resharp(tfs,mask.*R,voxelSize,ker_rad,tik_reg);
    mask_final = mask_ero;

    % save matrix
    if save_mat
        mkdir([path_mat '/rmbkg']);
        save([path_mat '/rmbkg/lfs_resharp.mat'],'lfs','-v7.3');
        save([path_mat '/mask/mask_resharp_final.mat'],'mask_final','-v7.3');
    end

    % save nifti
    mkdir([path_nft '/rmbkg/']);
    nii = make_nii(lfs,voxelSize);
    save_nii(nii,[path_nft '/rmbkg/lfs_resharp_xy.nii']);
    nii = make_nii(permute(lfs,[1 3 2]),voxelSize);
    save_nii(nii,[path_nft '/rmbkg/lfs_resharp_xz.nii']);
    nii = make_nii(permute(lfs,[2 3 1]),voxelSize);
    save_nii(nii,[path_nft '/rmbkg/lfs_resharp_yz.nii']);
    nii = make_nii(mask_final,voxelSize);
    save_nii(nii,[path_nft '/mask/mask_resharp_final.nii']);


    % inversion of susceptibility 
    disp('--> (9/9) TV susceptibility inversion on RE-SHARP...');
    sus = tvdi(lfs, mask_final, voxelSize, tv_reg, mag_cmb(:,:,:,4)); 
   
    % save matrix
    if save_mat
        mkdir([path_mat '/inversion']);
        save([path_mat '/inversion/sus_resharp.mat'],'sus','-v7.3');
    end

    % save nifti
    mkdir([path_nft '/inversion']);
    nii = make_nii(sus,voxelSize);
    save_nii(nii,[path_nft '/inversion/sus_resharp_xy.nii']);
    nii = make_nii(permute(sus,[1 3 2]),voxelSize);
    save_nii(nii,[path_nft '/inversion/sus_resharp_xz.nii']);
    nii = make_nii(permute(sus,[2 3 1]),voxelSize);
    save_nii(nii,[path_nft '/inversion/sus_resharp_yz.nii']);
end


%% E-SHARP
% if sum(strcmpi('esharp',bkgrm))
%     disp('--> (8/9) E-SHARP to remove background field ...');
%     Options.voxelSize = voxelSize;
%     lfs = esharp(tfs,mask.*R,Options);
%     mask_final = mask.*R;
% 
%     % save matrix
%     if save_mat
%         mkdir([path_mat '/rmbkg']);
%         save([path_mat '/rmbkg/lfs_esharp.mat'],'lfs','-v7.3');
%         save([path_mat '/mask/mask_esharp_final.mat'],'mask_final','-v7.3');
%     end
% 
%     % save nifti
%     mkdir([path_nft '/rmbkg/']);
%     nii = make_nii(lfs,voxelSize);
%     save_nii(nii,[path_nft '/rmbkg/lfs_esharp_xy.nii']);
%     nii = make_nii(permute(lfs,[1 3 2]),voxelSize);
%     save_nii(nii,[path_nft '/rmbkg/lfs_esharp_xz.nii']);
%     nii = make_nii(permute(lfs,[2 3 1]),voxelSize);
%     save_nii(nii,[path_nft '/rmbkg/lfs_esharp_yz.nii']);
%     nii = make_nii(mask_final,voxelSize);
%     save_nii(nii,[path_nft '/mask/mask_esharp_final.nii']);
% 
% 
%     % inversion of susceptibility 
%     disp('--> (9/9) TV susceptibility inversion on E-SHARP...');
%     sus = tvdi(lfs, mask_final, voxelSize, tv_reg, mag_cmb(:,:,:,4)); 
% 
%     % save matrix
%     if save_mat
%         mkdir([path_mat '/inversion']);
%         save([path_mat '/inversion/sus_esharp.mat'],'sus','-v7.3');
%     end
% 
%     % save nifti
%     mkdir([path_nft '/inversion']);
%     nii = make_nii(sus,voxelSize);
%     save_nii(nii,[path_nft '/inversion/sus_esharp_xy.nii']);
%     nii = make_nii(permute(sus,[1 3 2]),voxelSize);
%     save_nii(nii,[path_nft '/inversion/sus_esharp_xz.nii']);
%     nii = make_nii(permute(sus,[2 3 1]),voxelSize);
%     save_nii(nii,[path_nft '/inversion/sus_esharp_yz.nii']);
% end


%% PDF
if sum(strcmpi('pdf',bkgrm))
    disp('--> (8/9) PDF to remove background field ...');
    [lfs,mask_ero] = pdf(tfs,mask.*R,voxelSize,ker_rad,mag_cmb(:,:,:,4));

    % save matrix
    if save_mat
        mkdir([path_mat '/rmbkg']);
        save([path_mat '/rmbkg/lfs_pdf.mat'],'lfs','-v7.3');
        save([path_mat '/mask/mask_pdf_final.mat'],'mask_final','-v7.3');
    end

    % save nifti
    mkdir([path_nft '/rmbkg/']);
    nii = make_nii(lfs,voxelSize);
    save_nii(nii,[path_nft '/rmbkg/lfs_pdf_xy.nii']);
    nii = make_nii(permute(lfs,[1 3 2]),voxelSize);
    save_nii(nii,[path_nft '/rmbkg/lfs_pdf_xz.nii']);
    nii = make_nii(permute(lfs,[2 3 1]),voxelSize);
    save_nii(nii,[path_nft '/rmbkg/lfs_pdf_yz.nii']);
    nii = make_nii(mask_final,voxelSize);
    save_nii(nii,[path_nft '/mask/mask_pdf_final.nii']);


    % inversion of susceptibility 
    disp('--> (9/9) TV susceptibility inversion on PDF...');
    sus = tvdi(lfs, mask_final, voxelSize, tv_reg, mag_cmb(:,:,:,4)); 

    % save matrix
    if save_mat
        mkdir([path_mat '/inversion']);
        save([path_mat '/inversion/sus_pdf.mat'],'sus','-v7.3');
    end

    % save nifti
    mkdir([path_nft '/inversion']);
    nii = make_nii(sus,voxelSize);
    save_nii(nii,[path_nft '/inversion/sus_pdf_xy.nii']);
    nii = make_nii(permute(sus,[1 3 2]),voxelSize);
    save_nii(nii,[path_nft '/inversion/sus_pdf_xz.nii']);
    nii = make_nii(permute(sus,[2 3 1]),voxelSize);
    save_nii(nii,[path_nft '/inversion/sus_pdf_yz.nii']);
end
