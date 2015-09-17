function [lfs, mask_ero] = resharp(tfs,mask,vox,ker_rad,tik_reg, max_iter)
%RESHARP Background field removal
% Method is described in the following paper:
% Sun, H. and Wilman, A. H. (2013), 
% Background field removal using spherical mean value filtering and Tikhonov regularization. 
% Magn Reson Med. doi: 10.1002/mrm.24765
%
%   [LFS,MASK_ERO] = RESHARP(TFS,MASK,VOX,KER_RAD,TIK_REG,MAX_ITER)
%
%   LFS         : local field shift after background removal
%   MASK_ERO    : eroded mask after convolution
%   TFS         : input total field shift
%   MASK        : binary mask defining the brain ROI
%   VOX         : voxel size (mm), e.g. [1,1,1] for isotropic
%   KER_RAD     : radius of convolution kernel (mm), e.g. 5
%                 if length(ker_rad) ==2, a hollow spherical shell is used as kernel, 
%                 and these 2 parameters are the inner and outer radii
%   TIK_REG     : Tikhonov regularization parameter, e.g. 1e-3
%   MAX_ITER    : max # iterations for conjugate gradient solver (default: 200) 

if nargin < 6
    max_iter = 200 ;
end

imsize = size(tfs);

% make spherical/ellipsoidal convolution kernel (ker)
rx = round(ker_rad(1)/vox(1));
ry = round(ker_rad(1)/vox(2));
rz = round(ker_rad(1)/vox(3));
[X,Y,Z] = ndgrid(-rx:rx,-ry:ry,-rz:rz);

kernel = (X.^2)/(rx^2) + (Y.^2)/(ry^2) + (Z.^2)/(rz^2) <= 1 ;
nPointsKernel = sum(kernel(:));
kernel = kernel/nPointsKernel ;

% circularshift, linear conv to Fourier multiplication
csh = [rx,ry,rz]; % circularshift

% erode the mask by convolving with the kernel
cvsize = imsize + [2*rx+1, 2*ry+1, 2*rz+1] -1; % linear conv size
mask_tmp = real(ifftn(fftn(mask,cvsize).*fftn(kernel,cvsize)));
mask_tmp = mask_tmp(rx+1:end-rx, ry+1:end-ry, rz+1:end-rz); % same size
mask_ero = zeros(imsize);
mask_ero(mask_tmp > 1-1/nPointsKernel) = 1; % NO error tolerance

if length(ker_rad) == 2
    
    outerRadius = max(ker_rad) ;
    innerRadius = min(ker_rad) ;
    
    rxOut = round(outerRadius/vox(1));
    ryOut = round(outerRadius/vox(2));
    rzOut = round(outerRadius/vox(3));
    rxIn  = round(innerRadius/vox(1));
    ryIn  = round(innerRadius/vox(2));
    rzIn  = round(innerRadius/vox(3));
    
    kernelOut = (X.^2)/(rxOut^2) + (Y.^2)/(ryOut^2) + (Z.^2)/(rzOut^2) <= 1 ;
    kernelIn  = (X.^2)/(rxIn^2) + (Y.^2)/(ryIn^2) + (Z.^2)/(rzIn^2) <= 1 ;

    kernel    = kernelOut - kernelIn ;
    nPointsKernel = sum(kernel(:));

    kernel = kernel/nPointsKernel ;
end

% prepare convolution kernel: delta-kernel
dkernel = -kernel;
dkernel(rx+1,ry+1,rz+1) = 1-kernel(rx+1,ry+1,rz+1);
DKER = fftn(dkernel,imsize); % dkernel in Fourier domain


% RESHARP with Tikhonov regularization:   
%   argmin ||MSfCFx - MSfCFy||2 + lambda||x||2 
%   x: local field
%	y: total field
% 	M: binary mask
%	S: circular shift
%	F: forward FFT, f: inverse FFT (ifft) 
%	C: deconvolution kernelnel
%	lambda: tikhonov regularization parameter
%	||...||2: sum of square
%
%   create 'MSfCF' as an object 'H', then simplified as: 
%   argmin ||Hx - Hy||2 + lambda||x||2
%   To solve it, derivative equals 0: 
%   (H'H + lambda)x = H'Hy
%   Model as Ax = b, solve with cgs

H = cls_smvconv(imsize,DKER,csh,mask_ero); 
b = H'*(H*tfs(:));
m = cgs(@Afun, b, 1e-10, max_iter);

lfs = real(reshape(m,imsize)).*mask_ero;

% nested function
function y = Afun(x)
    y = H'*(H*x) + tik_reg*x;
    y = y(:);
end


end

