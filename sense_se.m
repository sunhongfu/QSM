function img_cmb = sense_se(img,reso,cref,radi)

% sense_se: SENSE combination for Single Echo
% Walsh paper to estimate coil sensitivities

% img: raw complex 4D images, 3D+receivers
% reso: resolution of the images (vector)
% cref: coil referece (coil number)
% radi: radius of the kernel (default: 3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reso = [0.5,0.752,2]; % (SWI 4.7T)
% cref = 3; % (3rd channel as reference coil for 4.7T)
% radi = 3; % (mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if no brain mask available (BET mask)
% temporarily SoS combine magnitude and extract ROI
mag_sos =  sum(abs(img).^2,4);
nii = make_nii(mag_sos,reso);
save_nii(nii,'mag_sos.nii');
% BET threshold relatively small to keep as much brain as possible  
! bet mag_sos.nii BET_sos -f 0.01 -m -R;  
! gunzip -f BET_sos.nii.gz BET_sos_mask.nii.gz;
nii = load_nii('BET_sos_mask.nii');
mask_sos = double(nii.img);

% image size: readout, phase encoding, slice encoding, receivers
[np,nv,nv2,nrcvrs] = size(img);


%% eigen value decomposition to estimate sensitivities
% define kernel (odd numbers)
kx = round(radi./reso(1));
ky = round(radi./reso(2));
kz = round(radi./reso(3));
kern = ones([kx*2+1,ky*2+1,kz*2+1]);

% construct coils correlation matrix
for j = 1:nrcvrs
    for k = 1:nrcvrs
        R(:,:,:,j,k) = img(:,:,:,j).*conj(img(:,:,:,k));
    end
end

% convolve with kernel (sum up in kernel)
cvsize = [np,nv,nv2] + [kx*2+1,ky*2+1,kz*2+1] -1;  % linear conv size
RS = ifft(ifft(ifft(  fft(fft(fft(R,cvsize(1),1),cvsize(2),2),cvsize(3),3) .* repmat(fftn(kern,cvsize), [1,1,1,nrcvrs,nrcvrs])  ,[],1),[],2),[],3);

clear R

% cut to the same size as before conv
RS = RS(1+kx:np+kx,  1+ky:nv+ky,  1+kz:nv2+kz, :,:);

RS = reshape(permute(RS,[4 5 1 2 3]),nrcvrs,nrcvrs,[]);

% extract only ROI of RS for eigen decomposition
RS_ex = RS(:,:,find(mask_sos(:)));

clear RS

sen_ex = zeros(nrcvrs,size(RS_ex,3));

matlabpool open
parfor i = 1:size(RS_ex,3)
    [V,D] = eig(RS_ex(:,:,i));
    sen_ex(:,i) = V(:,1);
end
matlabpool close

sen = zeros(nrcvrs,np*nv*nv2);
sen(:,find(mask_sos(:))) = sen_ex;
sen = reshape(permute(sen,[2,1]),[np,nv,nv2,nrcvrs]);

% relative to the ref coil
sen = sen./repmat(sen(:,:,:,cref)./abs(sen(:,:,:,cref)),[1 1 1 nrcvrs]);

clear RS_ex


%% combine coils using SENSE
img = reshape(img,[],nrcvrs);
sen = reshape(sen,[],nrcvrs);
img_cmb = sum(conj(sen).*img,2)./sum(sen.*conj(sen),2);
img_cmb = reshape(img_cmb,[np,nv,nv2]);
img_cmb(isnan(img_cmb)) = 0;
img_cmb(isinf(img_cmb)) = 0;


%% debug purpose
% phase after correction from each coil
nii = make_nii(angle(reshape(img,[np,nv,nv2,nrcvrs])./(reshape(sen,[np,nv,nv2,nrcvrs]))),reso);
save_nii(nii,'ph_cor.nii');

% magnitude and phase after combining
nii = make_nii(abs(img_cmb));
save_nii(nii,'mag_cmb.nii');
nii = make_nii(angle(img_cmb));
save_nii(nii,'ph_cmb.nii');
