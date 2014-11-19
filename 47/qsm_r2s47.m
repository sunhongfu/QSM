function qsm_r2s47(path_in, path_out, options)
%QSM_R2S47 Quantitative susceptibility mapping from R2* sequence at 4.7T.
%   QSM_R2S47(PATH_IN, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH_IN    - directory of .fid from gemsme3d sequence  : gemsme3d_R2s_01.fid
%   PATH_OUT   - directory to save nifti and/or matrixes   : QSM_R2s_vxxx
%   OPTIONS    - parameter structure including fields below
%    .bet_thr  - threshold for BET brain mask              : 0.5
%    .bkgrm    - background field removal method(s)        : 'resharp'
%    .smv_rad  - radius (mm) of SMV convolution kernel     : 6
%    .tik_reg  - Tikhonov regularization for RESHARP       : 0.001
%    .tsvd     - truncation of SVD for SHARP               : 0.05
%    .tv_reg   - Total variation regularization parameter  : 0.0005
%    .tvdi_n   - iteration number of TVDI (nlcg)           : 200
%    .echo_t   - keep only the first 'echo_t' echoes       : 5
%    .fit_t    - truncation level on fitting residual      : 5
%    .sav_all  - save all the variables for debug (~ 0)    : 1

%% default settings
if ~ exist('path_in','var') || isempty(path_in)
    path_in = pwd;
end

if exist([path_in '/fid'],'file')
    path_fid = path_in;
elseif exist([path_in '/gemsme3d_R2s_01.fid/fid'],'file')
    path_fid = [path_in, '/gemsme3d_R2s_01.fid'];
else
    error('cannot find .fid file');
end

if ~ exist('path_out','var') || isempty(path_out)
    path_out = path_fid;
end

if ~ exist('options','var') || isempty(options)
    options = [];
end

if ~ isfield(options,'bet_thr')
    options.bet_thr = 0.5;
end

if ~ isfield(options,'bkgrm')
    options.bkgrm = {'resharp','lbv'};
    % options.bkgrm = 'resharp';
end

if ~ isfield(options,'smv_rad')
    options.smv_rad = 6;
end

if ~ isfield(options,'tik_reg')
    options.tik_reg = 5e-4;
end

if ~ isfield(options,'tsvd')
    options.tsvd = 0.05;
end

if ~ isfield(options,'tv_reg')
    options.tv_reg = 5e-4;
end

if ~ isfield(options,'tvdi_n')
    options.tvdi_n = 200;
end

if ~ isfield(options,'echo_t')
    options.echo_t = 5;
end

if ~ isfield(options,'fit_t')
    options.fit_t = 10;
end

if ~ isfield(options,'sav_all')
    options.sav_all = 1;
end

bet_thr = options.bet_thr;
bkgrm   = options.bkgrm;
smv_rad = options.smv_rad;
tik_reg = options.tik_reg;
tsvd    = options.tsvd;
tv_reg  = options.tv_reg;
tvdi_n  = options.tvdi_n;
echo_t  = options.echo_t;
fit_t   = options.fit_t;
sav_all = options.sav_all;

%% define directories
path_qsm = [path_out '/QSM_R2s_v200_midnight'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);


%% reconstruct complex image from fid file
disp('--> reconstruct fid to complex img ...');
[img,par] = r2s47_recon(path_fid);


%%%%%%%%%%%%% old %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% keep only the first 'echo_t' echoes
% img_orig = img;
img = img(:,:,:,1:echo_t,:);
par.ne = echo_t;
TE = par.te + (0:par.ne-1)*par.esp;
%%%%%%%%%%%%% old %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% swap the first two dimensions and match scanner frame
img = permute(img,[2 1 3 4 5]); %PE,RO,SL,NE,RX
img = flipdim(img,1);
img = flipdim(img,2);
par.vox = [par.vox(2), par.vox(1), par.vox(3)];
[nv,np,nv2,ne,~] = size(img);
voxelSize = par.vox; % resolution in mm/pixel



phi = par.phi/180*pi;
psi = par.psi/180*pi;
theta = -par.theta/180*pi;
z_prjs = [sin(psi)*sin(theta), cos(psi)*sin(theta), cos(theta)]; 



%%%%%%%%%%%%%%%%%% new %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if par.nrcvrs > 1
    % combine using eig method
    img_cmb = coils_cmb(permute(img,[1 2 3 5 4]),voxelSize,TE,3,3);
else
    img_cmb = img;
end

mag_cmb = abs(img_cmb);
ph_cmb = angle(img_cmb);

mkdir('combine');
for echo = 1:size(img_cmb,4)
    nii = make_nii(mag_cmb(:,:,:,echo),voxelSize);
    save_nii(nii,['combine/mag_cmb' num2str(echo) '.nii']);
    nii = make_nii(ph_cmb(:,:,:,echo),voxelSize);
    save_nii(nii,['combine/ph_cmb' num2str(echo) '.nii']);
end


% %% generate mask from combined magnitude of the first echo
% disp('--> extract brain volume and generate mask ...');
% setenv('bet_thr',num2str(bet_thr));
% unix('bet combine/mag_cmb1.nii BET -f ${bet_thr} -m -R');
% unix('gunzip -f BET.nii.gz');
% unix('gunzip -f BET_mask.nii.gz');
% nii = load_nii('BET_mask.nii');
% mask = double(nii.img);


% % fit the echoes complexly
% [p1, dp1, relres, p0]=Fit_ppm_complex(img_cmb);
% nii = make_nii(p1,voxelSize);
% save_nii(nii,'complex_fit_p1.nii');


% % unwrap using laplacian
% unph = unwrapLaplacian(p1, size(p1), voxelSize);
% nii = make_nii(unph, voxelSize);
% save_nii(nii,'unph_lap.nii');

% %convert x to ppm
% tfs = unph/(2*pi*par.esp*42.576*4.7*1e6)*1e6;
% nii = make_nii(tfs, voxelSize);
% save_nii(nii,'tfs_ppm.nii');


% % RESHARP
% relres_blur = smooth3(relres,'box',round(6./par.res/2)*2+1); 
% nii = make_nii(relres_blur,voxelSize);
% save_nii(nii,'relres_blur.nii');

% R = ones(size(relres_blur));
% R(relres_blur >= 0.2) = 0;

% [lfs_resharp, mask_resharp] = resharp(tfs,mask.*R,voxelSize,smv_rad,tik_reg);
% mkdir('RESHARP');
% nii = make_nii(lfs_resharp,voxelSize);
% save_nii(nii,'RESHARP/lfs_resharp.nii');


% % inversion of susceptibility 
% disp('--> TV susceptibility inversion on RESHARP...');
% [sus_resharp,residual_resharp] = tvdi(lfs_resharp,mask_resharp,voxelSize,tv_reg,abs(img_cmb(:,:,:,end)),z_prjs,tvdi_n);


% % save nifti
% nii = make_nii(sus_resharp,voxelSize);
% save_nii(nii,'RESHARP/sus_resharp.nii');

% % save nifti
% nii = make_nii(sus_resharp.*mask_resharp,voxelSize);
% save_nii(nii,'RESHARP/sus_resharp_clean.nii');






%%%%%%%%%%%%%%%%%%%% old %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% combine magnitude channels using marc's 'arrayrec.m'
% disp('--> combine rcvrs for magnitude ...');
% mag_cmb = zeros(nv,np,nv2,ne);
% for echo = 1:ne
%     mag_cmb(:,:,:,echo) = arrayrec(squeeze(img(:,:,:,echo,:)),1/2);
% end


%% combine using eig
if par.nrcvrs > 1
    % combine using eig method
    mag_cmb = abs( sense_se(permute(img,[1 2 3 5 4]),voxelSize,3,3) );
else
    mag_cmb = abs(img);
end



% % save nifti for unwrapping usage
% mkdir('combine');
% for echo =  1:ne
%     nii = make_nii(mag_cmb(:,:,:,echo),voxelSize);
%     save_nii(nii,['combine/mag_te' num2str(echo) '.nii']);
% end


% % save nifti for unwrapping usage
% mkdir('combine');
% for echo =  1:ne
%     nii = make_nii(mag_cmb(:,:,:,echo),voxelSize);
%     save_nii(nii,['combine/mag_te' num2str(echo) '.nii']);
% end


%% generate mask from combined magnitude of the first echo
if ~ exist('BET_mask.nii','file')
    disp('--> extract brain volume and generate mask ...');
    setenv('bet_thr',num2str(bet_thr));
    unix('bet combine/mag_te1.nii BET -f ${bet_thr} -m -R');
    unix('gunzip -f BET.nii.gz');
    unix('gunzip -f BET_mask.nii.gz');
end
nii = load_nii('BET_mask.nii');
mask = double(nii.img);


% %% combine phase channels
% disp('--> combine rcvrs for phase ...');
% ph_cmb = sense_me(img,par);

% % save nifti for unwrapping usage
% for echo =  1:ne
%     nii = make_nii(ph_cmb(:,:,:,echo),voxelSize);
%     save_nii(nii,['combine/ph_te' num2str(echo) '.nii']);
% end

clear img;


%% unwrap phase from each echo
disp('--> unwrap aliasing phase for all TEs ...');

bash_command = sprintf(['for ph in combine/ph*\n' ...
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

unph_cmb = zeros(nv,np,nv2,ne);
for echo = 1:ne
    nii = load_nii(['unph_cmb' num2str(echo) '.nii']);
    unph_cmb(:,:,:,echo) = double(nii.img);
end


% check and correct for 2pi jump between echoes
disp('--> correct for potential 2pi jumps between TEs ...')

nii = load_nii('unph1.nii');
unph1 = double(nii.img);
nii = load_nii('unph2.nii');
unph2 = double(nii.img);
unph_diff = unph2 - unph1;

for echo = 2:ne
    meandiff = unph_cmb(:,:,:,echo)-unph_cmb(:,:,:,1)-(echo-1)*unph_diff;
    meandiff = meandiff(mask==1);
    meandiff = mean(meandiff(:));
    njump = round(meandiff/(2*pi));
    disp(['    ' num2str(njump) ' 2pi jumps for TE' num2str(echo)]);
    unph_cmb(:,:,:,echo) = unph_cmb(:,:,:,echo) - njump*2*pi;
    unph_cmb(:,:,:,echo) = unph_cmb(:,:,:,echo).*mask;
end


%% fit phase images with echo times
disp('--> magnitude weighted LS fit of phase to TE ...');
[tfs, fit_residual] = echofit(unph_cmb,mag_cmb,TE); 

% normalize to main field
% ph = gamma*dB*TE
% dB/B = ph/(gamma*TE*B0)
% units: TE s, gamma 2.675e8 rad/(sT), B0 4.7T
tfs = -tfs/(2.675e8*4.7)*1e6; % unit ppm


% generate reliability map
fit_residual_blur = smooth3(fit_residual,'box',round(6./par.res/2)*2+1); 
nii = make_nii(fit_residual_blur,voxelSize);
save_nii(nii,'fit_residual_blur.nii');

R = ones(size(fit_residual_blur));
R(fit_residual_blur >= 10) = 0;

% R = ones(size(fit_residual));
% R(fit_residual >= fit_t) = 0;


%% PDF
if sum(strcmpi('pdf',bkgrm))
    disp('--> PDF to remove background field ...');
    [lfs_pdf,mask_pdf] = pdf(tfs,mask.*R,voxelSize,smv_rad,mag_cmb(:,:,:,echo_t));

    % save nifti
    mkdir('PDF');
    nii = make_nii(lfs_pdf,voxelSize);
    save_nii(nii,'PDF/lfs_pdf.nii');

    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on PDF...');
    sus_pdf = tvdi(lfs_pdf, mask_pdf, voxelSize, tv_reg, mag_cmb(:,:,:,echo_t), z_prjs, tvdi_n); 

    % save nifti
    nii = make_nii(sus_pdf.*mask_pdf,voxelSize);
    save_nii(nii,'PDF/sus_pdf.nii');
end


%% SHARP (tsvd: truncation threthold for TSVD)
if sum(strcmpi('sharp',bkgrm))
    disp('--> SHARP to remove background field ...');
    [lfs_sharp, mask_sharp] = sharp(tfs,mask.*R,voxelSize,smv_rad,tsvd);

    % save nifti
    mkdir('SHARP');
    nii = make_nii(lfs_sharp,voxelSize);
    save_nii(nii,'SHARP/lfs_sharp.nii');
    
    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on SHARP...');
    sus_sharp = tvdi(lfs_sharp, mask_sharp, voxelSize, tv_reg, mag_cmb(:,:,:,echo_t), z_prjs, tvdi_n); 
   
    % save nifti
    nii = make_nii(sus_shar.*mask_sharp,voxelSize);
    save_nii(nii,'SHARP/sus_sharp.nii');
end


%% RE-SHARP (tik_reg: Tikhonov regularization parameter)
if sum(strcmpi('resharp',bkgrm))
    disp('--> RESHARP to remove background field ...');
    [lfs_resharp, mask_resharp] = resharp(tfs,mask.*R,voxelSize,smv_rad,tik_reg);


    % save nifti
    mkdir('RESHARP');
    nii = make_nii(lfs_resharp,voxelSize);
    save_nii(nii,'RESHARP/lfs_resharp.nii');


    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on RESHARP...');
    sus_resharp = tvdi(lfs_resharp, mask_resharp, voxelSize, tv_reg, mag_cmb(:,:,:,echo_t), z_prjs, tvdi_n); 
   

    % save nifti
    nii = make_nii(sus_resharp.*mask_resharp,voxelSize);
    save_nii(nii,'RESHARP/sus_resharp.nii');

end


%% LBV
if sum(strcmpi('lbv',bkgrm))
    disp('--> LBV to remove background field ...');
    lfs_lbv = LBV(tfs,mask.*R,size(tfs),voxelSize,0.01,2); % strip 2 layers
    mask_lbv = ones(size(mask));
    mask_lbv(lfs_lbv==0) = 0;

    % save nifti
    mkdir('lbv');
    nii = make_nii(lfs_lbv,voxelSize);
    save_nii(nii,'lbv/lfs_lbv.nii');


    % inversion of susceptibility 
    disp('--> TV susceptibility inversion on lbv...');
    sus_lbv = tvdi(lfs_lbv,mask_lbv,voxelSize,tv_reg,mag_cmb(:,:,:,echo_t),z_prjs,tvdi_n);   

    % save nifti
    nii = make_nii(sus_lbv.*mask_lbv,voxelSize);
    save_nii(nii,'lbv/sus_lbv.nii');

end

%% E-SHARP
% if sum(strcmpi('esharp',bkgrm))
%     disp('--> (8/9) E-SHARP to remove background field ...');
%     Options.voxelSize = voxelSize;
%     lfs = esharp(tfs,mask.*R,Options);
%     mask_final = mask.*R;
% 
% 
%     % save nifti
%     mkdir([path_qsm '/rmbkg/']);
%     nii = make_nii(lfs,voxelSize);
%     save_nii(nii,[path_qsm '/rmbkg/lfs_esharp_xy.nii']);
%     nii = make_nii(permute(lfs,[1 3 2]),voxelSize);
%     save_nii(nii,[path_qsm '/rmbkg/lfs_esharp_xz.nii']);
%     nii = make_nii(permute(lfs,[2 3 1]),voxelSize);
%     save_nii(nii,[path_qsm '/rmbkg/lfs_esharp_yz.nii']);
%     nii = make_nii(mask_final,voxelSize);
%     save_nii(nii,[path_qsm '/mask/mask_esharp_final.nii']);
% 
% 
%     % inversion of susceptibility 
%     disp('--> (9/9) TV susceptibility inversion on E-SHARP...');
%     sus = tvdi(lfs, mask_final, voxelSize, tv_reg, mag_cmb(:,:,:,4)); 
% 
% 
%     % save nifti
%     mkdir([path_qsm '/inversion']);
%     nii = make_nii(sus,voxelSize);
%     save_nii(nii,[path_qsm '/inversion/sus_esharp_xy.nii']);
%     nii = make_nii(permute(sus,[1 3 2]),voxelSize);
%     save_nii(nii,[path_qsm '/inversion/sus_esharp_xz.nii']);
%     nii = make_nii(permute(sus,[2 3 1]),voxelSize);
%     save_nii(nii,[path_qsm '/inversion/sus_esharp_yz.nii']);
% end
%%%%%%%%%%%%%%%%%%%% old %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% save all variables for debugging purpose
if sav_all
    clear nii;
    save('all.mat','-v7.3');
end

% save parameters used in the recon
save('parameters.mat','options','-v7.3')


%% clean up
% unix('rm *.nii*');
cd(init_dir);

