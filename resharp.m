function [lfs, mask_ero] = resharp(tfs,mask,par,radius,tik_reg)
%RESHARP Background field removal using RESHARP method.
%   [LSF,MASK_ERO] = RESHARP(TSF,MASK,PAR,TIK_REG)
%
%   TFS         : input total field shift
%   MASK        : binary mask defining the ROI
%   PAR         : parameter sets of the sequence
%   RADIUS      : Radius of convolution kernel size (mm)
%   TIK_REG     : Tikhonov regularization parameter
%   LFS         : local field shift after background removal
%   MASK_ERO    : eroded mask after convolution

imsize = size(tfs);
% resolution pix/mm
res = [par.lro/(par.np/2), par.lpe/par.nv, par.lpe2/par.nv2]*10; 

% make spherical/ellipsoidal convolution kernel
rx = round(radius/res(1));
ry = round(radius/res(2));
rz = round(radius/res(3));
[X,Y,Z] = ndgrid(-rx:rx,-ry:ry,-rz:rz);
h = (X.^2/rx^2 + Y.^2/ry^2 + Z.^2/rz^2 < 1);
ker = h/sum(h(:));

% circularshift, linear conv to Fourier multiplication
csh = [rx,ry,rz];

% erode the mask
mask_tmp = circshift(real(ifftn(fftn(mask).*fftn(ker,imsize))),-csh);
mask_ero = zeros(imsize);
mask_ero(mask_tmp > 1-6/sum(h(:))) = 1; % 6-1=5 error points tolerent


% prepare convolution kernel: delta-ker
dker = -ker;
dker(rx+1,ry+1,rz+1) = 1-ker(rx+1,ry+1,rz+1);
DKER = fftn(dker,imsize); % dker in Fourier domain


% RESHARP with Tikhonov regularization   
%   argmin ||MSfCFx - MSfCFy||2 + lambda||x||2 
%   x:local field, y:total field
%   create 'MSfCF' as an object 'H', then simplified as: 
%   argmin ||Hx - Hy||2 + lambda||x||2
%   to solve it, derivative equals 0: 
%   (H'H + lambda)x = H'Hy
%   Ax = b, solve with cgs

H = cls_smvconv(imsize,DKER,csh,mask_ero); 
b = H'*(H*tfs(:));
m = cgs(@Afun, b, 1e-6, 200);

lfs = real(reshape(m,imsize)).*mask_ero;


    function y = Afun(x)
        y = H'*(H*x) + tik_reg*x;
        y = y(:);
    end


end



