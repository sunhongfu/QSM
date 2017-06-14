function [lfs, mask_ero] = resharp_lsqr(tfs,mask,vox,ker_rad,iter_num)
%   [LSF,MASK_ERO,RES_TERM,REG_TERM] = RESHARP(TFS,MASK,VOX,KER_RAD,ITER_NUM)
%
%   LFS         : local field shift after background removal
%   MASK_ERO    : eroded mask after convolution
%   RES_TERM    : norm of data fidelity term
%   REG_TERM    : norm of regularization term
%   TFS         : input total field shift
%   MASK        : binary mask defining the brain ROI
%   VOX         : voxel size (mm), e.g. [1,1,1] for isotropic
%   KER_RAD     : radius of convolution kernel (mm), e.g. 3
%   ITER_NUM    : maximum number of CGS times of iteration, e.g. 200
%
%Method is described in the paper:
%Sun, H. and Wilman, A. H. (2013), 
%Background field removal using spherical mean value filtering and Tikhonov regularization. 
%Magn Reson Med. doi: 10.1002/mrm.24765

if ~ exist('vox','var') || isempty(vox)
    vox = [1,1,1];
end

if ~ exist('ker_rad','var') || isempty(ker_rad)
    ker_rad = 3;
end

if ~ exist('iter_num','var') || isempty(iter_num)
    iter_num = 500;
end


% make spherical/ellipsoidal convolution kernel (ker)
rx = round(ker_rad/vox(1));
ry = round(ker_rad/vox(2));
rz = round(ker_rad/vox(3));
rx = max(rx,2);
ry = max(ry,2);
rz = max(rz,2);
% rz = ceil(ker_rad/vox(3));
[X,Y,Z] = ndgrid(-rx:rx,-ry:ry,-rz:rz);
h = (X.^2/rx^2 + Y.^2/ry^2 + Z.^2/rz^2 <= 1);
ker = h/sum(h(:));

% circularshift, linear conv to Fourier multiplication
csh = [rx,ry,rz]; % circularshift

% pad zeros around to avoid errors in trans between linear conv and FT multiplication
tfs = padarray(tfs,csh);
mask = padarray(mask,csh);

imsize = size(tfs);

mask_ero = zeros(imsize);
mask_tmp = convn(mask,ker,'same');
%mask_ero(mask_tmp > 1-1/sum(h(:))) = 1; % no error points tolerence 
mask_ero(mask_tmp > 0.999999) = 1; % no error points tolerence 


% prepare convolution kernel: delta-ker
dker = -ker;
dker(rx+1,ry+1,rz+1) = 1-ker(rx+1,ry+1,rz+1);
DKER = fftn(dker,imsize); % dker in Fourier domain


% RESHARP with LSQR:   
%   MSfCFx = MSfCFy
%   x: local field
%	y: total field
% 	M: binary mask
%	S: circular shift
%	F: forward FFT, f: inverse FFT (ifft) 
%	C: deconvolution kernel

%   Model as Ax = b, solve with LSQR

b = mask_ero.*circshift(ifftn(DKER.*fftn(tfs)),-csh);
b = b(:);

m = lsqr(@Afun, b, 1e-6, iter_num);

lfs = real(reshape(m,imsize)).*mask_ero;

% remove the padding for outputs
lfs = lfs(rx+1:end-rx,ry+1:end-ry,rz+1:end-rz);
mask_ero = mask_ero(rx+1:end-rx,ry+1:end-ry,rz+1:end-rz);

% nested function

function y = Afun(x,transp_flag)
    if strcmp(transp_flag,'transp')
    % y = H'*(H*x) + tik_reg*x;
    x = reshape(x,imsize);
    y = ifftn(conj(DKER).*fftn(circshift(mask_ero.*x,csh)));

    y = y(:);
    
    else
    % y = H'*(H*x) + tik_reg*x;
    x = reshape(x,imsize);
    y = mask_ero.*circshift(ifftn(DKER.*fftn(x)),-csh);

    y = y(:); 
    end
end

end
