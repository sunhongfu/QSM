function [sus,residual_field] = tvdi(lfs, mask, vox, tv_reg, weights, z_prjs, Itnlim, pNorm)
%TVDI Total variation dipole inversion.

% Method is similar to Appendix in the following paper
% Lustig, M., Donoho, D. and Pauly, J. M. (2007), 
% Sparse MRI: The application of compressed sensing for rapid MR imaging. 
% Magn Reson Med, 58: 1182–1195. doi: 10.1002/mrm.21391

%   SUS = TVDI(LFS,MASK,VOX,TV_REG,WEIGHTS,THETA)
%
%   SUS    : susceptibility distribution after dipole inversion
%   LFS    : local field shift (field perturbation map)
%   MASK   : binary mask defining ROI
%   VOX    : voxel size, e.g. [1 1 1] for isotropic resolution
%   TV_REG : Total Variation regularization paramter, e.g. 5e-4
%   WEIGHTS: weights for the data consistancy term
%   Z_PRJS : normal vector of the imaging plane
%   ITNLIM : interation numbers of nlcg

if ~ exist('z_prjs','var') || isempty(z_prjs)
    z_prjs = [0, 0, 1]; % PURE axial slices
end

if ~ exist('Itnlim','var') || isempty(Itnlim)
    Itnlim = 500;
end

if ~ exist('pNorm','var') || isempty(pNorm)
    pNorm = 1;
end

[Nx,Ny,Nz] = size(lfs);
imsize = size(lfs);

% weights for data consistancy term (normalized)
W = mask.*weights;
W = W/sum(W(:))*sum(mask(:));

% % set the DC point of field in k-space to 0
% % mean value of lfs to be 0
% lfs = lfs.*mask;
% lfs = lfs-sum(lfs(:))/sum(mask(:));
% lfs = lfs.*mask;

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


% parameter structures for inversion
% data consistancy and TV term objects
param.FT = cls_dipconv([Nx,Ny,Nz],D);
% param.FT = cls_dipconv_mask([Nx,Ny,Nz],D,mask);
% param.FT = cls_dipconv_new([Nx,Ny,Nz],D,R);
param.TV = cls_tv;
% param.TV = cls_tv_mask(mask);

param.Itnlim = Itnlim; % interations numbers (adjust accordingly!)
param.gradToll = 1e-4; % step size tolerance stopping criterea
param.l1Smooth = eps; %1e-15; smoothing parameter of L1 norm
param.pNorm = pNorm; % type of norm to use (i.e. L1 L2 etc)
param.lineSearchItnlim = 100;
param.lineSearchAlpha = 0.01;
param.lineSearchBeta = 0.6;
param.lineSearchT0 = 1 ; % step size to start with

param.TVWeight = tv_reg; % TV penalty 
param.mask = mask; %%% not used in nlcg
param.data = lfs;

param.wt = W; % weighting matrix


tmp = fftn(lfs)./D;
T = 0.2; % truncation level
tmp(abs(D)<T) = 0;
sus_dc = real(mask.*ifftn(tmp));
nii = make_nii(sus_dc, vox);
save_nii(nii,'sus_dc.nii');

% non-linear conjugate gradient method
sus = nlcg(zeros(Nx,Ny,Nz), param);
% sus = nlcg(sus_dc, param);

% sus = real(sus).*mask;
% if want to keep the dipole fitting result
% don't mask it, instead, use the following:
sus = real(sus);

residual_field = lfs - real(ifftn(D.*fftn(sus)));

% nii = make_nii(real(ifftn(D.*fftn(sus))));
% save_nii(nii,'fit.nii');

% nii = make_nii(lfs,vox);
% save_nii(nii,'lfs.nii');
end
