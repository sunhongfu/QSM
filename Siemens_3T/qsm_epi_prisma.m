function qsm_epi_prisma(path_mag, path_ph, path_out, options)
%QSM_EPI_PRISMA Quantitative susceptibility mapping from EPI sequence at PRISMA (3T).
%   QSM_EPI_PRISMA(PATH_MAG, PATH_PH, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH_MAG     - directory of magnitude dicoms (mosaic)
%   PATH_PH      - directory of unfiltered phase dicoms (mosaic)
%   PATH_OUT     - directory to save nifti and/or matrixes   : QSM_EPI_PRISMA
%   OPTIONS      - parameter structure including fields below
%    .bet_thr    - threshold for BET brain mask              : 0.5
%    .bet_smooth - smoothness of BET brain mask at edges     : 2
%    .ph_unwrap  - 'prelude' or 'laplacian' or 'bestpath'    : 'prelude'
%    .bkg_rm     - background field removal method(s)        : 'resharp'
%                  options: 'pdf','sharp','resharp','esharp','lbv'
%                  to try all e.g.: {'pdf','sharp','resharp','esharp','lbv'}
%    .t_svd      - truncation of SVD for SHARP               : 0.1
%    .smv_rad    - radius (mm) of SMV convolution kernel     : 3
%    .tik_reg    - Tikhonov regularization for resharp       : 1e-3
%    .cgs_num    - max interation number for RESHARP         : 500
%    .lbv_peel   - LBV layers to be peeled off               : 2
%    .lbv_tol    - LBV interation error tolerance            : 0.01
%    .tv_reg     - Total variation regularization parameter  : 5e-4
%    .tvdi_n     - iteration number of TVDI (nlcg)           : 500

if ~ exist('path_mag','var') || isempty(path_mag)
    error('Please input the directory of magnitude DICOMs')
end

if ~ exist('path_ph','var') || isempty(path_ph)
    error('Please input the directory of unfiltered phase DICOMs')
end

if ~ exist('path_out','var') || isempty(path_out)
    path_out = pwd;
    display('Current directory for output')
end

if ~ exist('options','var') || isempty(options)
    options = [];
end

if ~ isfield(options,'bet_thr')
    options.bet_thr = 0.5;
end

if ~ isfield(options,'bet_smooth')
    options.bet_smooth = 2;
end

if ~ isfield(options,'ph_unwrap')
    options.ph_unwrap = 'prelude';
end

if ~ isfield(options,'bkg_rm')
    options.bkg_rm = 'resharp';
end

if ~ isfield(options,'t_svd')
    options.t_svd = 0.1;
end

if ~ isfield(options,'smv_rad')
    options.smv_rad = 3;
end

if ~ isfield(options,'tik_reg')
    options.tik_reg = 1e-3;
end

if ~ isfield(options,'cgs_num')
    options.cgs_num = 500;
end

if ~ isfield(options,'lbv_tol')
    options.lbv_tol = 0.01;
end

if ~ isfield(options,'lbv_peel')
    options.lbv_peel = 2;
end

if ~ isfield(options,'tv_reg')
    options.tv_reg = 5e-4;
end

if ~ isfield(options,'inv_num')
    options.inv_num = 500;
end


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


% read in DICOMs of both magnitude and raw unfiltered phase images
path_mag = cd(cd(path_mag));
path_ph = cd(cd(path_ph));
path_out = cd(cd(path_out));
mag_list = dir(path_mag);
ph_list = dir(path_ph);
mag_list = mag_list(~strncmpi('.', {mag_list.name}, 1));
ph_list = ph_list(~strncmpi('.', {ph_list.name}, 1));

% mosaic form to 4D nifti
for i = 1:length(mag_list)
    mag_mosaic(:,:,i) = single(dicomread([path_mag,filesep,mag_list(i).name]));
end
for i = 1:length(ph_list)
    ph_mosaic(:,:,i) = single(dicomread([path_ph,filesep,ph_list(i).name]));
    ph_mosaic(:,:,i) = ph_mosaic(:,:,i)/4095*2*pi - pi;
end

% crop mosaic into individual images
dicom_info = dicominfo([path_mag,filesep,mag_list(1).name]);
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

mag_all = zeros(wRow,wCol,nSL,size(mag_mosaic,3));
ph_all = mag_all;
for i = 1:size(mag_mosaic,3)
    for x = 1:wRow
        for y = 1:wCol
            for z = 1:nSL
                X = floor((z-1)/nCol)*wRow + x;
                Y = mod(z-1,nCol)*wCol + y;
                mag_all(x,y,z,i) = mag_mosaic(X,Y,i);
                ph_all(x,y,z,i) = ph_mosaic(X,Y,i);
            end
        end
    end
end

% permute the images to 
% x:right-to-left
% y:anterior-to-posterior
% z:inferior-to-superior
mag_all = permute(mag_all,[2 1 3 4]);
ph_all = permute(ph_all,[2 1 3 4]);

% get the sequence parameters
vox = [dicom_info.PixelSpacing(1), dicom_info.PixelSpacing(2), dicom_info.SliceThickness];
imsize = size(mag_all);
[~,~,~,nVol] = size(mag_all);

% angles!!! (z projections)
Xz = dicom_info.ImageOrientationPatient(3);
Yz = dicom_info.ImageOrientationPatient(6);
Zxyz = cross(dicom_info.ImageOrientationPatient(1:3),dicom_info.ImageOrientationPatient(4:6));
Zz = Zxyz(3);
%Zz = sqrt(1 - Xz^2 - Yz^2);
z_prjs = [Xz, Yz, Zz];


% define output directories
path_qsm = [path_out '/QSM_EPI_PRISMA'];
mkdir(path_qsm);
init_dir = pwd;
cd(path_qsm);

% save magnitude and phase for each echo
mkdir('src');
for i = 1:size(mag_all,4) % all time series
    nii = make_nii(mag_all(:,:,:,i),vox);
    save_nii(nii,['src/mag' num2str(i,'%03i') '.nii']);
    nii = make_nii(ph_all(:,:,:,i),vox);
    save_nii(nii,['src/ph' num2str(i,'%03i') '.nii']);
end


% generate BET on 1st volume
[status,cmdout] = unix('rm BET*');
disp('--> extract brain volume and generate mask ...');
setenv('bet_thr',num2str(bet_thr));
setenv('bet_smooth',num2str(bet_smooth));
bash_script = ['bet2 src/mag001.nii BET ' ...
    '-f ${bet_thr} -m -w ${bet_smooth}'];
unix(bash_script);
unix('gunzip -f BET.nii.gz');
unix('gunzip -f BET_mask.nii.gz');
nii = load_nii(['BET_mask.nii']);
mask = double(nii.img);


if nVol > 1
    % generate masks for all volumes
    mask_rep = repmat(mask,[1 1 1 size(mag_all,4)]);
    nii = make_nii(mask_rep,vox);
    save_nii(nii,'mask_rep.nii');

    nii = make_nii(mag_all,vox);
    save_nii(nii,'mag_all.nii');
    
    % spm to align all volumes
    P = spm_select('ExtList', pwd, '^mag_all.nii',Inf);
    flags.mask=0;
    spm_realign(P);
    load mag_all.mat
    m=[vox(1) 0 0 0; 0 vox(2) 0 0; 0 0 vox(3) 0; 0 0 0 1];
    for i = 2:size(mat,3)
        % mat(:,:,i) = inv(inv(m)*mat(:,:,i))*m;
            mat(:,:,i) = m\(mat(:,:,i))*m;
    end
    save('mask_rep.mat','mat');
    P = spm_select('ExtList', pwd, '^mask_rep.nii',Inf);
    flags.mask=0;
    spm_reslice(P,flags);

    nii = load_nii('rmask_rep.nii');
    mask_all = nii.img;
    mask_all(isnan(mask_all)) = 0;
    mask_all(isinf(mask_all)) = 0;
else
    mask_all = mask;
end

% do the remaining QSM steps individually for each volume
for i = 1:nVol % all time series
    nii = make_nii(mask_all(:,:,:,i),vox);
    save_nii(nii,['BET' num2str(i,'%03i') '_mask.nii']);
    mask = mask_all(:,:,:,i);
    mag = mag_all(:,:,:,i);
    ph = ph_all(:,:,:,i);
    
    % unwrap the phase
    if strcmpi('prelude',ph_unwrap)
        % unwrap phase with PRELUDE
        disp('--> unwrap aliasing phase ...');
        setenv('time_series',num2str(i,'%03i'));
        bash_script = ['prelude -a src/mag${time_series}.nii ' ...
            '-p src/ph${time_series}.nii -u unph${time_series}.nii ' ...
            '-m BET${time_series}_mask.nii -n 8'];
        unix(bash_script);
        unix('gunzip -f unph${time_series}.nii.gz');
        nii = load_nii(['unph' num2str(i,'%03i') '.nii']);
        unph = double(nii.img);

    elseif strcmpi('laplacian',ph_unwrap)
        % Ryan Topfer's Laplacian unwrapping
        disp('--> unwrap aliasing phase using laplacian...');
        Options.voxelSize = vox;
        unph = lapunwrap(ph, Options);
        nii = make_nii(unph, vox);
        save_nii(nii,['unph_lap' num2str(i,'%03i') '.nii']);

    elseif strcmpi('bestpath',ph_unwrap)
        % unwrap the phase using best path
        disp('--> unwrap aliasing phase using bestpath...');
            [pathstr, ~, ~] = fileparts(which('3DSRNCP.m'));
        setenv('pathstr',pathstr);
        setenv('nv',num2str(imsize(1)));
        setenv('np',num2str(imsize(2)));
        setenv('ns',num2str(imsize(3)));

        fid = fopen('wrapped_phase.dat','w');
        fwrite(fid,ph,'float');
        fclose(fid);
        % mask_unwrp = uint8(hemo_mask.*255);
        mask_unwrp = uint8(abs(mask)*255);
        fid = fopen('mask_unwrp.dat','w');
        fwrite(fid,mask_unwrp,'uchar');
        fclose(fid);

        bash_script = ['${pathstr}/3DSRNCP wrapped_phase.dat mask_unwrp.dat unwrapped_phase.dat ' ...
            '$nv $np $ns reliability.dat'];
        unix(bash_script);

        fid = fopen('unwrapped_phase.dat','r');
        tmp = fread(fid,'float');
        unph = reshape(tmp - round(mean(tmp(mask==1))/(2*pi))*2*pi,imsize(1:3)).*mask;
        fclose(fid);

        nii = make_nii(unph, vox);
        save_nii(nii,['unph_bestpath' num2str(i,'%03i') '.nii']);
    end

	% normalize to ppm unit
	tfs = unph/(2.675e8*dicom_info.EchoTime*dicom_info.MagneticFieldStrength)*1e9; % unit ppm


    % background field removal
    % PDF
    if sum(strcmpi('pdf',bkg_rm))
        disp('--> PDF to remove background field ...');
        [lfs_pdf,mask_pdf] = projectionontodipolefields(tfs,mask,vox,smv_rad,mag,z_prjs);
        % 2D 2nd order polyfit to remove any residual background
        lfs_pdf = lfs_pdf - poly2d(lfs_pdf,mask_pdf);

        % save nifti
        mkdir('PDF');
        nii = make_nii(lfs_pdf,vox);
        save_nii(nii,['PDF/lfs_pdf' num2str(i,'%03i') '.nii']);

        % inversion of susceptibility 
        disp('--> TV susceptibility inversion on PDF...');
        sus_pdf = tvdi(lfs_pdf,mask_pdf,vox,tv_reg,mag,z_prjs,inv_num); 

        % demean susceptibility
        sus_pdf = sus_pdf - mean(sus_pdf(logical(mask_pdf(:))));

        % save nifti
        nii = make_nii(sus_pdf,vox);
        save_nii(nii,['PDF/sus_pdf' num2str(i,'%03i') '.nii']);
    end


    % SHARP (t_svd: truncation threthold for t_svd)
    if sum(strcmpi('sharp',bkg_rm))
        disp('--> SHARP to remove background field ...');
        [lfs_sharp, mask_sharp] = sharp(tfs,mask,vox,smv_rad,t_svd);
        % 2D 2nd order polyfit to remove any residual background
        lfs_sharp = lfs_sharp - poly2d(lfs_sharp,mask_sharp);

        % save nifti
        mkdir('SHARP');
        nii = make_nii(lfs_sharp,vox);
        save_nii(nii,['SHARP/lfs_sharp' num2str(i,'%03i') '.nii']);
        
        % inversion of susceptibility 
        disp('--> TV susceptibility inversion on SHARP...');
        sus_sharp = tvdi(lfs_sharp,mask_sharp,vox,tv_reg,mag,z_prjs,inv_num); 

        % demean susceptibility
        sus_sharp = sus_sharp - mean(sus_sharp(logical(mask_sharp(:))));

        % save nifti
        nii = make_nii(sus_sharp.*mask_sharp,vox);
        save_nii(nii,['SHARP/sus_sharp' num2str(i,'%03i') '.nii']);
    end


    % RE-SHARP (tik_reg: Tikhonov regularization parameter)
    if sum(strcmpi('resharp',bkg_rm))
    	disp('--> RESHARP to remove background field ...');
        [lfs_resharp, mask_resharp] = resharp(tfs,mask,vox,smv_rad,tik_reg,cgs_num);
        % 2D 2nd order polyfit to remove any residual background
        lfs_resharp = lfs_resharp - poly2d(lfs_resharp,mask_resharp);

        % save nifti
        mkdir('RESHARP');
        nii = make_nii(lfs_resharp,vox);
        save_nii(nii,['RESHARP/lfs_resharp' num2str(i,'%03i') '.nii']);

        % inversion of susceptibility 
        disp('--> TV susceptibility inversion on RESHARP...');
        sus_resharp = tvdi(lfs_resharp,mask_resharp,vox,tv_reg,mag,z_prjs,inv_num); 
        
        % demean susceptibility
        sus_resharp = sus_resharp - mean(sus_resharp(logical(mask_resharp(:))));

        % save nifti
        nii = make_nii(sus_resharp.*mask_resharp,vox);
        save_nii(nii,['RESHARP/sus_resharp' num2str(i,'%03i') '.nii']);
    end


    % E-SHARP (SHARP edge extension)
    if sum(strcmpi('esharp',bkg_rm))
        disp('--> E-SHARP to remove background field ...');
        Parameters.voxelSize             = vox; % in mm
        Parameters.resharpRegularization = tik_reg ;
        Parameters.resharpKernelRadius   = smv_rad ; % in mm
        Parameters.radius                = [ 10 10 5 ] ;

        % pad matrix size to even number
        pad_size = mod(size(tfs),2);
        tfs = padarray(tfs.*mask, pad_size, 'post');

        % taking off additional 1 voxels from edge - not sure the outermost 
        % phase data included in the original mask is reliable. 
        mask_shaved = shaver( ( tfs ~= 0 ), 1 ) ; % 1 voxel taken off
        totalField  = mask_shaved .* tfs ;

        % resharp 
        [reducedLocalField, maskReduced] = ...
            resharp( totalField, ...
                     double(mask_shaved), ...
                     Parameters.voxelSize, ...
                     Parameters.resharpKernelRadius, ...
                     Parameters.resharpRegularization ) ;

        % extrapolation ~ esharp 
        reducedBackgroundField = maskReduced .* ( totalField - reducedLocalField) ;

        extendedBackgroundField = extendharmonicfield( ...
           reducedBackgroundField, mask, maskReduced, Parameters) ;

        backgroundField = extendedBackgroundField + reducedBackgroundField ;
        localField      = totalField - backgroundField ;

        lfs_esharp      = localField(1+pad_size(1):end,1+pad_size(2):end,1+pad_size(3):end);
        mask_esharp     = mask_shaved(1+pad_size(1):end,1+pad_size(2):end,1+pad_size(3):end);  

        % 2D 2nd order polyfit to remove any residual background
        lfs_esharp = lfs_esharp - poly2d(lfs_esharp,mask_esharp);

        % save nifti
        mkdir('ESHARP');
        nii = make_nii(lfs_esharp,vox);
        save_nii(nii,['ESHARP/lfs_esharp' num2str(i,'%03i') '.nii']);

        % inversion of susceptibility 
        disp('--> TV susceptibility inversion on ESHARP...');
        sus_esharp = tvdi(lfs_esharp,mask_esharp,vox,tv_reg,mag,z_prjs,inv_num); 

       % demean susceptibility
        sus_esharp = sus_esharp - mean(sus_esharp(logical(mask_esharp(:))));

        % save nifti
        nii = make_nii(sus_esharp.*mask_esharp,vox);
        save_nii(nii,['ESHARP/sus_esharp' num2str(i,'%03i') '.nii']);
    end

    % LBV
    if sum(strcmpi('lbv',bkg_rm))
        disp('--> LBV to remove background field ...');
        lfs_lbv = LBV(tfs,mask,imsize,vox,lbv_tol,lbv_peel); % strip 2 layers
        mask_lbv = ones(size(mask));
        mask_lbv(lfs_lbv==0) = 0;
        % 2D 2nd order polyfit to remove any residual background
        lfs_lbv = lfs_lbv - poly2d(lfs_lbv,mask_lbv);

        % save nifti
        mkdir('LBV');
        nii = make_nii(lfs_lbv,vox);
        save_nii(nii,['LBV/lfs_lbv' num2str(i,'%03i') '.nii']);

        % inversion of susceptibility 
        disp('--> TV susceptibility inversion on lbv...');
        sus_lbv = tvdi(lfs_lbv,mask_lbv,vox,tv_reg,mag,z_prjs,inv_num);   

       % demean susceptibility
        sus_lbv = sus_lbv - mean(sus_lbv(logical(mask_lbv(:))));

        % save nifti
        nii = make_nii(sus_lbv.*mask_lbv,vox);
        save_nii(nii,['LBV/sus_lbv' num2str(i,'%03i') '.nii']);
    end

end


if nVol > 1
    % align all susceptiblity maps using the matrix from magnitude
    % read in all susceptiblity nii into 4D
    if sum(strcmpi('pdf',bkg_rm))
        mkdir('PDF/spm_realign')
        listing=dir('PDF/sus_pdf*.nii');
        for j = 1:length(listing)
            nii = load_nii(['PDF',filesep,listing(j).name]);
            sus_pdf_all(:,:,:,j) = double(nii.img);
        end
        nii = make_nii(sus_pdf_all,vox);
        save_nii(nii,'PDF/spm_realign/sus_pdf_all.nii');

        % align sus_pdf_all
        unix('cp mag_all.mat PDF/spm_realign/sus_pdf_all.mat');
        cd('PDF/spm_realign');
        P = spm_select('ExtList', pwd, '^sus_pdf_all.nii',Inf);
        flags.mask=0;
        spm_reslice(P,flags);
        cd ../..
    end

    if sum(strcmpi('sharp',bkg_rm))
        mkdir('SHARP/spm_realign')
        listing=dir('SHARP/sus_sharp*.nii');
        for j = 1:length(listing)
            nii = load_nii(['SHARP',filesep,listing(j).name]);
            sus_sharp_all(:,:,:,j) = double(nii.img);
        end
        nii = make_nii(sus_sharp_all,vox);
        save_nii(nii,'SHARP/spm_realign/sus_sharp_all.nii');

        % align sus_sharp_all
        unix('cp mag_all.mat SHARP/spm_realign/sus_sharp_all.mat');
        cd('SHARP/spm_realign');
        P = spm_select('ExtList', pwd, '^sus_sharp_all.nii',Inf);
        flags.mask=0;
        spm_reslice(P,flags);
        cd ../..
    end

    if sum(strcmpi('resharp',bkg_rm))
        mkdir('RESHARP/spm_realign')
        listing=dir('RESHARP/sus_resharp*.nii');
        for j = 1:length(listing)
            nii = load_nii(['RESHARP',filesep,listing(j).name]);
            sus_resharp_all(:,:,:,j) = double(nii.img);
        end
        nii = make_nii(sus_resharp_all,vox);
        save_nii(nii,'RESHARP/spm_realign/sus_resharp_all.nii');

        % align sus_resharp_all
        unix('cp mag_all.mat RESHARP/spm_realign/sus_resharp_all.mat');
        cd('RESHARP/spm_realign');
        P = spm_select('ExtList', pwd, '^sus_resharp_all.nii',Inf);
        flags.mask=0;
        spm_reslice(P,flags);
        cd ../..
    end

    if sum(strcmpi('esharp',bkg_rm))
        mkdir('ESHARP/spm_realign')
        listing=dir('ESHARP/sus_esharp*.nii');
        for j = 1:length(listing)
            nii = load_nii(['ESHARP',filesep,listing(j).name]);
            sus_esharp_all(:,:,:,j) = double(nii.img);
        end
        nii = make_nii(sus_esharp_all,vox);
        save_nii(nii,'ESHARP/spm_realign/sus_esharp_all.nii');

        % align sus_esharp_all
        unix('cp mag_all.mat ESHARP/spm_realign/sus_esharp_all.mat');
        cd('ESHARP/spm_realign');
        P = spm_select('ExtList', pwd, '^sus_esharp_all.nii',Inf);
        flags.mask=0;
        spm_reslice(P,flags);
        cd ../..
    end

    if sum(strcmpi('lbv',bkg_rm))
        mkdir('LBV/spm_realign')
        listing=dir('LBV/sus_lbv*.nii');
        for j = 1:length(listing)
            nii = load_nii(['LBV',filesep,listing(j).name]);
            sus_lbv_all(:,:,:,j) = double(nii.img);
        end
        nii = make_nii(sus_lbv_all,vox);
        save_nii(nii,'LBV/spm_realign/sus_lbv_all.nii');

        % align sus_lbv_all
        unix('cp mag_all.mat LBV/spm_realign/sus_lbv_all.mat');
        cd('LBV/spm_realign');
        P = spm_select('ExtList', pwd, '^sus_lbv_all.nii',Inf);
        flags.mask=0;
        spm_reslice(P,flags);
        cd ../..
    end
end

save('all.mat','-v7.3');
cd(init_dir);
