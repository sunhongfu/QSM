function img_cmb_all = coils_cmb(img,vox,cref,radi,te,off_corr)

% D. Walsh paper to estimate coil sensitivities
% Adaptive reconstruction of phased array MR imagery. MRM 2000

% img: raw complex 4D or 5D images: [3Dimages_dimensions, multi-receivers_dimension, (multi-echoes_dimension)]
% vox: resolution of the images (vector), voxel size
% te: at least the first and second TEs are needed: [te(1), te(2)]
% cref: coil referece (coil number)
% radi: radius of the kernel (e.g. 3mm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for example:
% vox = [1 ,1, 1];  % isotropic 1mm resolution
% te = [3, 7, 11, 15, 19]; % echo times, only the first two are required!
%   if single-echo data, 'te' will not be used, just replace 'te' with [] as input
% cref = 1; % (1st channel as reference coil)
% radi = 4; % (mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ exist('cref','var') || isempty(cref)
    cref = 1;
end

if ~ exist('radi','var') || isempty(radi)
    radi = 4;
end

if ~ exist('te','var') || isempty(te)
    te = [];
end

if ~ exist('off_corr','var') || isempty(off_corr)
    off_corr = 0;
end


% image size: readout, phase encoding, slice encoding, receivers, echoes
[np,nv,nv2,nrcvrs,ne] = size(img);
img_orig = img;
img = img(:,:,:,:,1);


%% eigen value decomposition to estimate sensitivities
% define kernel (odd numbers)
kx = round(radi./vox(1));
ky = round(radi./vox(2));
kz = round(radi./vox(3));
kern = ones([kx*2+1,ky*2+1,kz*2+1]);

% construct coils correlation matrix
R = zeros([np,nv,nv2,nrcvrs,nrcvrs]);
for j = 1:nrcvrs
    for k = 1:nrcvrs
        R(:,:,:,j,k) = img(:,:,:,j).*conj(img(:,:,:,k));
    end
end

% convolve with kernel (sum up in kernel)
cvsize = [np,nv,nv2] + [kx*2+1,ky*2+1,kz*2+1] -1;  % linear conv size
RS = ifft(ifft(ifft(  fft(fft(fft(R,cvsize(1),1),cvsize(2),2),cvsize(3),3) ...
    .* repmat(fftn(kern,cvsize), [1,1,1,nrcvrs,nrcvrs])  ,[],1),[],2),[],3);

clear R

% cut to the same size as before conv
RS = RS(1+kx:np+kx,  1+ky:nv+ky,  1+kz:nv2+kz, :,:);
RS = reshape(permute(RS,[4 5 1 2 3]),nrcvrs,nrcvrs,np*nv,nv2);


sen = zeros(nrcvrs,np*nv,nv2);
matlabpool open
parfor sl = 1:nv2
    sen(:,:,sl) = eig_fun(RS(:,:,:,sl));
end
matlabpool close
sen = reshape(permute(sen,[2,3,1]),[np,nv,nv2,nrcvrs]);

% relative to the ref coil
sen = sen./repmat(sen(:,:,:,cref)./abs(sen(:,:,:,cref)),[1 1 1 nrcvrs]);

clear RS_ex


%% combine coils using SENSE
img_cmb_all = zeros([np,nv,nv2,ne]);
sen = reshape(sen,[],nrcvrs);
for echo = 1:ne
    img = reshape(img_orig(:,:,:,:,echo),[],nrcvrs);
    img_cmb = sum(conj(sen).*img,2)./sum(sen.*conj(sen),2);
    img_cmb = reshape(img_cmb,[np,nv,nv2]);
    img_cmb(isnan(img_cmb)) = 0;
    img_cmb(isinf(img_cmb)) = 0;
    img_cmb_all(:,:,:,echo) = img_cmb;
end



%% if multi-echo then correct for phase-offset
if (echo > 1) && off_corr
    nii = make_nii(abs(img_cmb_all(:,:,:,1)),vox);
    save_nii(nii,'mag1.nii');
    nii = make_nii(abs(img_cmb_all(:,:,:,2)),vox);
    save_nii(nii,'mag2.nii');
    nii = make_nii(angle(img_cmb_all(:,:,:,1)),vox);
    save_nii(nii,'wrph1.nii');
    nii = make_nii(angle(img_cmb_all(:,:,:,2)),vox);
    save_nii(nii,'wrph2.nii');

    %% generate mask from combined magnitude of the first echo
    disp('--> extract brain volume and generate mask ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bet_thr = 0.5; % adjust the BET parameter accordingly, refer to BET from FSL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    setenv('bet_thr',num2str(bet_thr));
    unix('bet mag1.nii BET -f ${bet_thr} -m -R');
    unix('gunzip -f BET.nii.gz');
    unix('gunzip -f BET_mask.nii.gz');

    bash_command = sprintf(['prelude -a mag1.nii -p wrph1.nii -u unph1.nii -m BET_mask -n 8&\n' ...
                            'prelude -a mag2.nii -p wrph2.nii -u unph2.nii -m BET_mask -n 8&\n' ...
                            'wait\n' ...
                            'gunzip -f unph*.gz\n']);
    unix(bash_command);
    nii = load_nii('unph1.nii');
    unph1 = double(nii.img);
    nii = load_nii('unph2.nii');
    unph2 = double(nii.img);

    offset = (te(1)*unph2 - te(2)*unph1)/(te(1)-te(2));
    img_cmb_all = img_cmb_all./exp(1j.*repmat(offset,[1,1,1,echo]));
end

end







function sen = eig_fun(RS)
    for i = 1:size(RS,3);
        [V,D] = eig(RS(:,:,i));
        sen(:,i) = V(:,1);
    end
end
