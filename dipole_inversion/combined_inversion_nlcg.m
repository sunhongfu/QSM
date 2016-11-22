function [chi, res] = combined_inversion_nlcg(tfs, mask, vox, z_prjs, ker_rad, Tik_weight, TV_weight)

% || m * F_{-1} * C * D * F * chi - m * F_{-1} * C * F * B_t|| 
% 	+ Tik_weight * ||m * chi||
%   + TV_weight * TV(m * chi)


[Nx,Ny,Nz] = size(tfs);

% create K-space filter kernel D
%%%%% make this a seperate function in the future
FOV = vox.*[Nx,Ny,Nz];
FOVx = FOV(1);
FOVy = FOV(2);
FOVz = FOV(3);

x = -Nx/2:Nx/2-1;
y = -Ny/2:Ny/2-1;
z = -Nz/2:Nz/2-1;
[kx,ky,kz] = ndgrid(x/FOVx,y/FOVy,z/FOVz);
% D = 1/3 - kz.^2./(kx.^2 + ky.^2 + kz.^2);
D = 1/3 - (kx.*z_prjs(1)+ky.*z_prjs(2)+kz.*z_prjs(3)).^2./(kx.^2 + ky.^2 + kz.^2);
D(floor(Nx/2+1),floor(Ny/2+1),floor(Nz/2+1)) = 0;
D = fftshift(D);



% REHARP part
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
mask_ero = zeros([Nx,Ny,Nz]);
mask_tmp = convn(mask,ker,'same');
mask_ero(mask_tmp > 1-1/sum(h(:))) = 1; % no error points tolerence 


% prepare convolution kernel: delta-ker
dker = -ker;
dker(rx+1,ry+1,rz+1) = 1-ker(rx+1,ry+1,rz+1);
DKER = fftn(dker,[Nx,Ny,Nz]); % dker in Fourier domain


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter structures for inversion
% data consistancy and TV term objects
params.FT = cls_singleStep(mask_ero,DKER,D,csh);
params.TV =	cls_tv;

params.Itnlim = 50; % interations numbers (adjust accordingly!)
params.gradToll = 1e-6; % step size tolerance stopping criterea
params.l1Smooth = eps; %1e-15; smoothing parameter of L1 norm
params.pNorm = 1; % type of norm to use (i.e. L1 L2 etc)
params.lineSearchItnlim = 100;
params.lineSearchAlpha = 0.01;
params.lineSearchBeta = 0.6;
params.lineSearchT0 = 1 ; % step size to start with

params.Tik_weight = Tik_weight; 
params.TV_weight = TV_weight; % TV penalty 
params.mask = mask_ero; %%% not used in nlcg
params.data = mask_ero.*circshift(ifftn(DKER.*fftn(tfs)),-csh);
% params.wt = mask.*mag(:,:,:,end); % weighting matrix
params.wt = mask_ero; % weighting matrix



% non-linear conjugate gradient method
chi = nlcg_singlestep(zeros(Nx,Ny,Nz), params);

% if want to keep the dipole fitting result
% don't mask it, instead, use the following:
% chi = real(chi).*mask;
chi = real(chi);

% residual difference between fowardly calculated field and lfs
res = tfs - real(ifftn(D.*fftn(chi)));


