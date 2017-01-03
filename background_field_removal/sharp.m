function [lfs, mask_ero] = sharp(tfs,mask,vox,ker_rad,tsvd)
%SHARP Background field removal
%   [LSF,MASK_ERO] = SHARP(TFS,MASK,VOX,KER_RAD,TSVD)
%
%   LFS         : local field shift after background removal
%   MASK_ERO    : eroded mask after convolution
%   TFS         : input total field shift
%   MASK        : binary mask defining the brain ROI
%   VOX         : voxel size (mm), e.g. [1,1,1] for isotropic
%   KER_RAD     : radius of convolution kernel (mm), e.g. 5
%   TSVD        : truncated singular value decomposition. e.g. 0.05

if ~ exist('vox','var') || isempty(vox)
    vox = [1,1,1];
end

if ~ exist('ker_rad','var') || isempty(ker_rad)
    ker_rad = 4;
end

if ~ exist('tik_reg','var') || isempty(tik_reg)
    tik_reg = 5e-4;
end

imsize = size(tfs);

% make spherical/ellipsoidal convolution kernel (ker)
rx = round(ker_rad/vox(1));
ry = round(ker_rad/vox(2));
rz = round(ker_rad/vox(3));
rx = max(rx,1);
ry = max(ry,1);
rz = max(rz,1);
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


% SHARP
% convolute the total field (ext + int) with d_kernel
ph_tmp = circshift(ifftn(fftn(tfs).*DKER),-csh);
% erode the result (abandon brain edges)
ph_tmp = ph_tmp.*mask_ero;
% deconvolution
ph_int = fftn(ph_tmp)./DKER;
ph_int(abs(DKER)<tsvd) = 0;
ph_tmp = circshift(ifftn(ph_int),csh);
lfs = real(ph_tmp).*mask_ero;

