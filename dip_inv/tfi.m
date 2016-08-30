function [chi,polyfit,res] = tfi(tfs, Res_wt, sus_mask, Tik_mask, TV_mask, Tik_reg, TV_reg, vox, z_prjs, Itnlim)

% argmin ||Res_wt * (F_{-1} * D * F * sus_mask * chi + V *p - tfs)|| + Tik_reg*||Tik_mask * chi|| + TV_reg*TV(TV_mask * chi)


% if ~ exist('weights','var') || isempty(weights)
%     weights = mask_b;
% end

% normalize the weights
Res_wt = Res_wt/sqrt(sum(Res_wt(:).^2)/numel(Res_wt));


[Nx,Ny,Nz] = size(tfs);
[np nv nv2] = size(tfs);

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


% create polynomials (2nd order 3d for now)
poly_order = 2;

px = repmat((1:np)',[nv*nv2,1]);
py = repmat((1:nv),[np,1]);
py = repmat(py(:),[nv2,1]);
pz = repmat((1:nv2),[np*nv,1]);
pz = pz(:);

if poly_order == 1
	% first order polyfit
	P = [px, py, pz, ones(length(px),1)];
	% P_nz = [px_nz, py_nz, pz_nz, ones(length(px_nz),1)]; % polynomials

elseif poly_order == 2	
	% second order
	P = [px.^2, py.^2, pz.^2, px.*py, px.*pz, py.*pz, px, py, pz, ones(length(px),1)]; % polynomials
	% P_nz = [px_nz.^2, py_nz.^2, pz_nz.^2, px_nz.*py_nz, px_nz.*pz_nz, py_nz.*pz_nz, px_nz, py_nz, pz_nz, ones(length(px_nz),1)]; % polynomials

elseif poly_order == 3
	% third order
	P = [px.^3, py.^3, pz.^3, px.*py.^2, px.*pz.^2, px.*py.*pz, py.*px.^2, py.*pz.^2, pz.*px.^2, pz.*py.^2 ...
		px.^2, py.^2, pz.^2, px.*py, px.*pz, py.*pz, px, py, pz, ones(length(px),1)]; % polynomials
	% P_nz = [px_nz.^3, py_nz.^3, pz_nz.^3, px_nz.*py_nz.^2, px_nz.*pz_nz.^2, px_nz.*py_nz.*pz_nz, py_nz.*px_nz.^2, py_nz.*pz_nz.^2, pz_nz.*px_nz.^2, pz_nz.*py_nz.^2 ...
		% px_nz.^2, py_nz.^2, pz_nz.^2, px_nz.*py_nz, px_nz.*pz_nz, py_nz.*pz_nz, px_nz, py_nz, pz_nz, ones(length(px_nz),1)]; % polynomials
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter structures for inversion
% data consistancy and TV term objects
params.FT = cls_dipconv_polyfit([Nx,Ny,Nz],D,sus_mask,P);
params.TV = cls_tv;

params.Itnlim = Itnlim; % interations numbers (adjust accordingly!)
params.gradToll = 1e-6; % step size tolerance stopping criterea
params.l1Smooth = eps; %1e-15; smoothing parameter of L1 norm
params.pNorm = 1; % type of norm to use (i.e. L1 L2 etc)
params.lineSearchItnlim = 100;
params.lineSearchAlpha = 0.01;
params.lineSearchBeta = 0.6;
params.lineSearchT0 = 1 ; % step size to start with

params.Tik_reg = Tik_reg; 
params.TV_reg = TV_reg; 
params.Tik_mask = Tik_mask; 
params.TV_mask = TV_mask; 
params.sus_mask = sus_mask;
params.Res_wt = Res_wt;
params.data = tfs;
params.imsize = [Nx,Ny,Nz];


% initial estimate
% (1) polyfit
polyfit = poly3d(tfs,Tik_mask,2);
coeff = P\polyfit(:);
D0 = D;
D0(abs(D)<0.25) = 0.25;
chi0 = real(ifftn(fftn((tfs-polyfit).*Tik_mask)./D0)).*Tik_mask;
nii = make_nii(chi0,vox);
save_nii(nii,'chi0.nii');

% non-linear conjugate gradient method
[chi, coeff, Res_term, TV_term, Tik_term] = nlcg_dipconv_polyfit([zeros(Nx*Ny*Nz,1);coeff], params);

% if want to keep the dipole fitting result
% don't mask it, instead, use the following:
% chi = real(chi).*mask;
chi = real(chi);

polyfit = reshape(P*coeff,[Nx,Ny,Nz]);

% residual difference between fowardly calculated field and lfs
res = tfs - real(ifftn(D.*fftn(chi)));

