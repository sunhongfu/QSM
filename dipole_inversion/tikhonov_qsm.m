function [chi, Res_term, TV_term, Tik_term] = tikhonov_qsm(tfs, Res_wt, sus_mask, TV_mask, Tik_mask, air_mask, TV_reg, Tik_reg, TV_reg2, vox, P, z_prjs, Itnlim)
% argmin ||Res_wt * (F_{-1} * D * F * sus_mask * chi - tfs)|| + TV_reg * TV|TV_mask * chi| + Tik_reg * ||Tik_mask * chi|| + TV_reg2 * TV|air_mask * chi|
%
% tfs:      total field shift; can be local field shift if use this for local field inversion
% Res_wt:   weighting matrix for the residual/fidelity term, usually brain mask
% sus_mask: input 1 for total field inversion; input brain mask for local field inversion
% TV_mask:  usually brain mask for TV regularization of the local tissue susceptiblity distribution
% Tik_mask: usually brain mask for Tikhonov regularization of the local tissue susceptiblity distribution
% TV_reg:   Regularization parameter for TV term
% Tik_reg:  Regularization parameter for Tikhonov term
% vox:      voxel size in mm
% z_prjs:   z-projections of the k-space coordinates onto main field B0 direction
% Itnlim:   Iteration number limit 

if ~ exist('P','var') || isempty(P)
    P = 1; % no preconditioning
end

if ~ exist('z_prjs','var') || isempty(z_prjs)
    z_prjs = [0, 0, 1]; % PURE axial slices
end

if ~ exist('Itnlim','var') || isempty(Itnlim)
    Itnlim = 500;
end


% normalize the weights
% Res_wt = Res_wt/sqrt(sum(Res_wt(:).^2)/numel(Res_wt));

% Res_wt = TV_mask.*Res_wt;
% Res_wt = Res_wt/sum(Res_wt(:))*sum(TV_mask(:));

% create K-space filter kernel D
%%%%% make this a seperate function in the future
[Nx, Ny, Nz] = size(tfs);
FOV  = vox.*[Nx, Ny, Nz];
FOVx = FOV(1);
FOVy = FOV(2);
FOVz = FOV(3);

x = -Nx/2:Nx/2-1;
y = -Ny/2:Ny/2-1;
z = -Nz/2:Nz/2-1;

[kx, ky, kz] = ndgrid(x/FOVx,y/FOVy,z/FOVz);
D = 1/3 - (kx.*z_prjs(1) + ky.*z_prjs(2) + kz.*z_prjs(3)).^2 ./ (kx.^2 + ky.^2 + kz.^2);
D(floor(Nx/2+1), floor(Ny/2+1), floor(Nz/2+1)) = 0;
D = fftshift(D);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter structures for inversion
% data consistancy and TV term objects
params.FT               = cls_dipconv([Nx, Ny, Nz], D); % class for dipole kernel convolution
params.TV               = cls_tv; 						% class for TV operation

params.Itnlim           = Itnlim; 						% interations numbers (adjust accordingly!)
params.gradToll         = 1e-4; 						% step size tolerance stopping criterea
params.l1Smooth         = eps; 							% 1e-15; smoothing parameter of L1 norm
params.pNorm            = 1; 							% type of norm to use (i.e. L1 L2 etc)
params.lineSearchItnlim = 100;
params.lineSearchAlpha  = 0.01;
params.lineSearchBeta   = 0.6;
params.lineSearchT0     = 1 ; 							% step size to start with

params.Tik_reg          = Tik_reg; 
params.TV_reg           = TV_reg; 
params.Tik_mask         = Tik_mask; 
params.TV_mask          = TV_mask; 
params.sus_mask         = sus_mask;
params.Res_wt           = Res_wt;
params.data             = tfs;
params.P                = P;

%params.P                = sus_mask.*(Tik_mask + 60*(1-Tik_mask));
%params.P                = sus_mask.*(Tik_mask + 60*(1-Tik_mask)) + (1-sus_mask);

% params.air_mask = 1;
% params.TV_reg2 = 0;
params.air_mask = air_mask;
params.TV_reg2 = TV_reg2;

% non-linear conjugate gradient method
[chi, Res_term, TV_term, Tik_term] = nlcg_tik(zeros([Nx, Ny, Nz]), params);


% LSQR method
% argmin ||Res_wt * (F_{-1} * D * F * sus_mask * chi - tfs)|| + Tik_reg * ||Tik_mask * chi|| 

% b = P.*sus_mask.*ifftn(D.*fftn(Res_wt.*Res_wt.*tfs));
% b = b(:);
% imsize = size(tfs);
% m = lsqr(@Afun, b, 1e-6, Itnlim);
% chi = real(reshape(m,imsize));

% function y = Afun(x,transp_flag)
%     if strcmp(transp_flag,'transp')
%     % y = H'*(H*x) + tik_reg*x;
%     x = reshape(x,imsize);
%     y = P.*sus_mask.*ifftn(D.*fftn(Res_wt.*Res_wt.*ifftn(D.*fftn(sus_mask.*P.*x)))) + Tik_reg*P.*Tik_mask.*Tik_mask.*P.*x;

%     y = y(:);
    
%     else
%     % y = H'*(H*x) + tik_reg*x;
%     x = reshape(x,imsize);
%     y = P.*sus_mask.*ifftn(D.*fftn(Res_wt.*Res_wt.*ifftn(D.*fftn(sus_mask.*P.*x)))) + Tik_reg*P.*Tik_mask.*Tik_mask.*P.*x;

%     y = y(:); 
%     end
% end


% if want to keep the dipole fitting result
% don't mask it, instead, use the following:
chi = real(params.P.*chi);
% otherwise mask it
end
