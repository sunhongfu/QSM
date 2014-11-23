function [ph_corr, mask_ero] = sharp(ph,mask,res,radius,tsvd)
%SHARP Background field removal using SHARP method.
%   [PH_SHP,MASK_ERO] = sharp(PH,MASK,PAR,TSVD)
%
%   PH:       input total field
%   MASK:     binary mask defining the ROI
%   RES:      resolution of the images (vector)
%   TSVD:     truncation level for TSVD
%   PH_SHP:   local field after background removal (SHARP)
%   MASK_ERO: new ROI after erosion

imsize = size(ph);

% make spherical/ellipsoidal convolution kernel
rx = round(radius/res(1));
ry = round(radius/res(2));
rz = round(radius/res(3));
[X,Y,Z] = ndgrid(-rx:rx,-ry:ry,-rz:rz);
h = (X.^2/rx^2 + Y.^2/ry^2 + Z.^2/rz^2 < 1);
ker = h/sum(h(:));

% circularshift, linear conv to Fourier multiplication
csh = [rx,ry,rz]; % circularshift

% erode the mask by convolving with the kernel
cvsize = imsize + [2*rx+1, 2*ry+1, 2*rz+1] -1; % linear conv size
mask_tmp = real(ifftn(fftn(mask,cvsize).*fftn(ker,cvsize)));
mask_tmp = mask_tmp(rx+1:end-rx, ry+1:end-ry, rz+1:end-rz); % same size
mask_ero = zeros(imsize);
mask_ero(mask_tmp > 1-6/sum(h(:))) = 1; % 5 points error tolerance

% prepare convolution kernel: delta-ker
dker = -ker;
dker(rx+1,ry+1,rz+1) = 1-ker(rx+1,ry+1,rz+1);
D_ker = fftn(dker,imsize); % dker in Fourier domain

% SHARP
% convolute the total field (ext + int) with d_kernel
ph_tmp = circshift(ifftn(fftn(ph).*D_ker),-csh);
% erode the result (abandon brain edges)
ph_tmp = ph_tmp.*mask_ero;
% deconvolution
PH_int = fftn(ph_tmp)./D_ker;
PH_int(abs(D_ker)<tsvd) = 0;
ph_tmp = circshift(ifftn(PH_int),csh);
ph_corr = real(ph_tmp).*mask_ero;

