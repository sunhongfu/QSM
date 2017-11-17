function [img_cmb_all,sen] = adaptive_cmb(img,vox,cref,radi,flag_pool)

% D. Walsh paper to estimate coil sensitivities
% Adaptive reconstruction of phased array MR imagery. MRM 2000

% img: raw complex 4D or 5D images: [3Dimages, multi-receivers, (multi-echoes)]
% vox: resolution of the images (vector), voxel size
% cref: coil referece (coil number)
% radi: radius of the kernel (e.g. 5mm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vox = [1 ,1, 1];  % isotropic 1mm resolution
% cref = 1; % (1st channel as reference coil)
% radi = 5; % (mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ exist('vox','var') || isempty(vox)
    vox = [1 1 1];
end

if ~ exist('cref','var') || isempty(cref)
    cref = 1;
end

if ~ exist('radi','var') || isempty(radi)
    radi = 5;
end

if ~ exist('flag_pool','var') || isempty(flag_pool)
    flag_pool = 1;
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
R = zeros([np,nv,nv2,nrcvrs,nrcvrs],'single');
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
% (1) using MATLAB POOL
if flag_pool
    if exist('parpool')
        poolobj=parpool;
    else
        matlabpool open
    end
    parfor sl = 1:nv2
        sen(:,:,sl) = eig_fun(RS(:,:,:,sl));
    end
    if exist('parpool')
        delete(poolobj);
    else
        matlabpool close
    end
else
% (2) NOT using MATLAB POOL
    for sl = 1:nv2
        sen(:,:,sl) = eig_fun(RS(:,:,:,sl));
    end
end

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

end







function sen = eig_fun(RS)
    for i = 1:size(RS,3);
        % [V,D] = eig(RS(:,:,i));
        [V,D] = svd(RS(:,:,i));
        sen(:,i) = V(:,1);
    end
end
