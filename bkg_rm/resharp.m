function [lfs, mask_ero, data_fidelity, regularization_term] = resharp(tfs,mask,vox,ker_rad,tik_reg,cgs_num)
%   [LSF,MASK_ERO] = RESHARP(TFS,MASK,VOX,KER_RAD,TIK_REG)
%
%   LFS         : local field shift after background removal
%   MASK_ERO    : eroded mask after convolution
%   TFS         : input total field shift
%   MASK        : binary mask defining the brain ROI
%   VOX         : voxel size (mm), e.g. [1,1,1] for isotropic
%   KER_RAD     : radius of convolution kernel (mm), e.g. 3
%   TIK_REG     : Tikhonov regularization parameter, e.g. 1e-3
%	CGS_NUM     : maximum number of CGS times of iteration
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

if ~ exist('tik_reg','var') || isempty(tik_reg)
    tik_reg = 1e-3;
end

if ~ exist('cgs_num','var') || isempty(cgs_num)
    cgs_num = 500;
end

imsize = size(tfs);

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

% erode the mask by convolving with the kernel
% cvsize = imsize + [2*rx+1, 2*ry+1, 2*rz+1] -1; % linear conv size
% mask_tmp = real(ifftn(fftn(mask,cvsize).*fftn(ker,cvsize)));
% mask_tmp = mask_tmp(rx+1:end-rx, ry+1:end-ry, rz+1:end-rz); % same size
mask_ero = zeros(imsize);
mask_tmp = convn(mask,ker,'same');
mask_ero(mask_tmp > 1-1/sum(h(:))) = 1; % no error points tolerence 


% prepare convolution kernel: delta-ker
dker = -ker;
dker(rx+1,ry+1,rz+1) = 1-ker(rx+1,ry+1,rz+1);
DKER = fftn(dker,imsize); % dker in Fourier domain


% RESHARP with Tikhonov regularization:   
%   argmin ||MSfCFx - MSfCFy||2 + lambda||x||2 
%   x: local field
%	y: total field
% 	M: binary mask
%	S: circular shift
%	F: forward FFT, f: inverse FFT (ifft) 
%	C: deconvolution kernel
%	lambda: tikhonov regularization parameter
%	||...||2: sum of square
%
%   create 'MSfCF' as an object 'H', then simplified as: 
%   argmin ||Hx - Hy||2 + lambda||x||2
%   To solve it, derivative equals 0: 
%   (H'H + lambda)x = H'Hy
%   Model as Ax = b, solve with cgs

% H = cls_smvconv(imsize,DKER,csh,mask_ero); 

b = ifftn(conj(DKER).*fftn(circshift(mask_ero.*circshift(ifftn(DKER.*fftn(tfs)),-csh),csh)));
b = b(:);

% b = H'*(H*tfs(:));
m = cgs(@Afun, b, 1e-6, cgs_num);

lfs = real(reshape(m,imsize)).*mask_ero;

data_fidelity = mask_ero.*circshift(ifftn(DKER.*fftn(tfs-reshape(m,imsize))),-csh);
data_fidelity = norm(data_fidelity(:));

regularization_term = norm(m);

% nested function
function y = Afun(x)
    % y = H'*(H*x) + tik_reg*x;
    x = reshape(x,imsize);
    y = ifftn(conj(DKER).*fftn(circshift(mask_ero.*circshift(ifftn(DKER.*fftn(x)),-csh),csh))) + tik_reg*x;

    y = y(:);
end


end

