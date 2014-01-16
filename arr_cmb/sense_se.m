function img_cmb_all = sense_se(img,vox,cref,radi)

% sense_se: SENSE combination for Single Echo 
% Walsh paper to estimate coil sensitivities

% img: raw complex 4D images, 3D+receivers (+echoes)
% vox: resolution of the images (vector)
% cref: coil referece (coil number)
% radi: radius of the kernel (default: 3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vox = [0.5,0.752,2]; % (SWI 4.7T)
% cref = 3; % (3rd channel as reference coil for 4.7T)
% radi = 3; % (mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% image size: readout, phase encoding, slice encoding, receivers
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

% nii = make_nii(abs(sen),vox);
% save_nii(nii,'sen_mag.nii');
% nii = make_nii(angle(sen),vox);
% save_nii(nii,'sen_ph.nii');

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
        [V,D] = eig(RS(:,:,i));
        sen(:,i) = V(:,1);
    end
end
