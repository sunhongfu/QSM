function [ph_corr, mask_ero] = sharp(ph,mask,par,radius,t)
%SHARP Background field removal using SHARP method.
%   [PH_SHP,MASK_ERO] = sharp(PH,MASK,PAR,TSVD)
%
%   PH:       input total field
%   MASK:     binary mask defining the ROI
%   PAR:      parameter sets of the sequence
%   TSVD:     truncation level for TSVD
%   PH_SHP:   local field after background removal (SHARP)
%   MASK_ERO: new ROI after erosion


% t: TSVD threshold

% ph can be either 3D (np*nv*ns) or 4D (np*nv*ns*ne) or 5D (np*nv*ns*ne*nrcvrs)
[np nv ns] = size(ph);
res = [par.lro/(par.np/2), par.lpe/par.nv, par.lpe2/par.nv2]*10; % resolution pix/mm

% make spherical/ellipsoidal convolution kernel
rx = round(radius/res(1));
ry = round(radius/res(2));
rz = round(radius/res(3));
[X,Y,Z] = ndgrid(-rx:rx,-ry:ry,-rz:rz);
h = (X.^2/rx^2 + Y.^2/ry^2 + Z.^2/rz^2 < 1);
ker = h/sum(h(:));

% erode the mask
mask_tmp = circshift(real(ifftn(fftn(mask).*fftn(ker,[np nv ns]))),[-rx,-ry,-rz]);
mask_ero = zeros([np nv ns]);
mask_ero(mask_tmp > 1-6/sum(h(:))) = 1; % 2 error points tolerance

% prepare deconvolution kernel
d_ker = -ker;
d_ker(rx+1,ry+1,rz+1) = 1-ker(rx+1,ry+1,rz+1);
D_ker = fftn(d_ker,[np nv ns]); % d_ker in Fourier domain

% SHARP
% convolute the total field (ext + int) with d_kernel
ph_tmp = circshift(ifftn(fftn(ph).*D_ker),[-rx,-ry,-rz]);
% erode the result (abandon brain edges)
ph_tmp = ph_tmp.*mask_ero;
% deconvolution
PH_int = fftn(ph_tmp)./D_ker;
PH_int(abs(D_ker)<t) = 0;
ph_tmp = circshift(ifftn(PH_int),[rx,ry,rz]);
ph_corr = real(ph_tmp).*mask_ero;


% % %%
% % % set the DC point of field in k-space to 0
% % k = fftn(ph_corr);
% % k(1,1,1) = 0;
% % ph_corr = real(ifftn(k)).*mask_ero;
% % %%
