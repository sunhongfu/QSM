function sus = tvdi(lfs, mask, par, TVWeight)
%TVDI Total variation dipole inversion.
%   SUS = tvdi(LFS,MASK,PAR,TV_REG)
%
%   LFS    : local field shift (field perturbation map)
%   MASK   : binary mask defining ROI
%   PAR    : parameter set for the sequence
%   TV_REG : Total variation regularization paramter
%   SUS    : susceptibility distribution after dipole inversion


[np,nv,ns] = size(lfs);

W1 = ones([np nv ns]); % weights for data consistancy term

% % normalize lfs
% % ph = gamma*dB*TE
% % dB/B = ph/(gamma*TE*B0)
% % units: TE s, gamma 2.675e8 rad/(sT), B0 4.7T
lfs = -lfs/(2.675e8*4.7)*1e6; % unit ppm

% set the DC point of field in k-space to 0
lfs = lfs-mean(lfs(:));

%%%%%%%%%%%%%%%%%%%%% needs scaling with FOV %%%%%%%%%%%%%%%%%%%%%%%
% create K-space filter kernel D
FOVx = par.lro;
FOVy = par.lpe;
FOVz = par.lpe2;

x = -np/2:np/2-1;
y = -nv/2:nv/2-1;
z = -ns/2:ns/2-1;
[kx,ky,kz] = ndgrid(x/FOVx,y/FOVy,z/FOVz);
D = 1/3 - kz.^2./(kx.^2 + ky.^2 + kz.^2);

D = fftshift(D);
D(1,1,1) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.FT = cls_dipconv([np,nv,ns],D);
param.TV = cls_tv;

param.TVWeight = TVWeight;     % TV penalty 
param.mask = mask;

param.Itnlim = 200;
param.gradToll = 0;     % step size tolerance stopping criterea
param.l1Smooth = eps; %1e-15; 	% smoothing parameter of L1 norm
param.pNorm = 1;            % type of norm to use (i.e. L1 L2 etc)
param.lineSearchItnlim = 100;
param.lineSearchAlpha = 0.01;
param.lineSearchBeta = 0.6;
param.lineSearchT0 = 1 ;    % step size to start with
param.data = lfs;
param.wt = W1; % weighting matrix

% tmp = fftn(lfs)./D;
% T = 0.1; % truncation level for initial guess
% tmp(abs(D)<T) = 0;
% sus_dc = real(mask.*ifftn(tmp));
% nii = make_nii(sus_dc);
% save_nii(nii,['sus_threshold_' num2str(T) '.nii']);


% non-linear conjugate gradient method
sus = nlcg(zeros(np,nv,ns), param);
sus = real(sus).*mask;


end
