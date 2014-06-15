function sus = tvdi(lfs, mask, vox, tv_reg, weights, nor_vec, Itnlim)
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
%   NOR_VEC: normal vector of the imaging plane
%   ITNLIM : interation numbers of nlcg

if ~ exist('nor_vec','var') || isempty(nor_vec)
    nor_vec = [0, 0, 1]; % PURE axial slices
end

if ~ exist('Itnlim','var') || isempty(Itnlim)
    Itnlim = 200;
end

[Nx,Ny,Nz] = size(lfs);

% weights for data consistancy term (normalized)
W = mask.*weights;
W = W/sum(W(:))*sum(mask(:));

% set the DC point of field in k-space to 0
% mean value of lfs to be 0
lfs = lfs.*mask;
lfs = lfs-sum(lfs(:))/sum(mask(:));
lfs = lfs.*mask;

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
D = 1/3 - kz.^2./(kx.^2 + ky.^2 + kz.^2);

D(floor(Nx/2+1),floor(Ny/2+1),floor(Nz/2+1)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rotate the kernel to match the angle of acqisition
% the following is just an example of how to rotate
% it has to be adjusted to rotate accordingly

% theta = -acos(nor_vec(3));

% if abs(theta) >= 0.01
%     disp('rotate dipole kernel along first dimension axis');
%     [X,Y,Z] = ndgrid(x,y,z);
%     rotM = [1,0,0; 0,cos(theta),-sin(theta); 0,sin(theta),cos(theta)];

%     afterR = rotM*[X(:)' ; Y(:)' ; Z(:)'];
%     finalM = afterR + repmat([Nx/2+1;Ny/2+1;Nz/2+1],[1 Nx*Ny*Nz]);
%     Index = finalM;

%     IndexX = Index(1,:);
%     IndexY = Index(2,:);
%     IndexZ = Index(3,:);

%     IndexX(IndexX<=1) = 1;
%     IndexY(IndexY<=1) = 1;
%     IndexZ(IndexZ<=1) = 1;

%     IndexX(IndexX>Nx) = Nx;
%     IndexY(IndexY>Ny) = Ny;
%     IndexZ(IndexZ>Nz) = Nz;

%     % bilinear interpolation
%     D_r = zeros(size(D));
%     for i = 1:length(IndexX)
%         ix = IndexX(i);
%         iy = IndexY(i);
%         iz = IndexZ(i);
%         d = iz - floor(iz);
%         d1 = iy - floor(iy);
%         D_r(i) = D(ix,floor(iy),floor(iz))*(1-d)*(1-d1) + D(ix,floor(iy),ceil(iz))*d*(1-d1)...
%             + D(ix,ceil(iy),floor(iz))*d1*(1-d) + D(ix,ceil(iy),ceil(iz))*d*d1;
%     end

%     D = D_r;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rotate the kernel
if ~ isequal(nor_vec,[0 0 1])
    R = rotationmatrix([0 0 1], nor_vec); % Rotation Matrix, from [0 0 1] to nor_vec
   % R = rotationmatrix(nor_vec,[0 0 1]);
	[X,Y,Z] = ndgrid(x,y,z);
	afterR = R*[X(:)' ; Y(:)' ; Z(:)']; % rotate about the 0 origin
	finalM = afterR + repmat([Nx/2+1;Ny/2+1;Nz/2+1],[1 Nx*Ny*Nz]); % translation
	Index = finalM;

	% new Index corresponding to old one
	IndexX = Index(1,:);
	IndexY = Index(2,:);
	IndexZ = Index(3,:);

	IndexX(IndexX<=1) = 1;
	IndexY(IndexY<=1) = 1;
	IndexZ(IndexZ<=1) = 1;

	IndexX(IndexX>Nx) = Nx;
	IndexY(IndexY>Ny) = Ny;
	IndexZ(IndexZ>Nz) = Nz;


	% D kernel after rotation
	D_r = interpn(D,IndexX,IndexY,IndexZ,'cubic');
	D_r = reshape(D_r,size(D));

	D = D_r;
	nii = make_nii(D,vox);
	save_nii(nii,'D.nii');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


D = fftshift(D);


% parameter structures for inversion
% data consistancy and TV term objects
param.FT = cls_dipconv([Nx,Ny,Nz],D);
param.TV = cls_tv;

param.Itnlim = Itnlim; % interations numbers (adjust accordingly!)
param.gradToll = 0; % step size tolerance stopping criterea
param.l1Smooth = eps; %1e-15; smoothing parameter of L1 norm
param.pNorm = 1; % type of norm to use (i.e. L1 L2 etc)
param.lineSearchItnlim = 100;
param.lineSearchAlpha = 0.01;
param.lineSearchBeta = 0.6;
param.lineSearchT0 = 1 ; % step size to start with

param.TVWeight = tv_reg; % TV penalty 
param.mask = mask;
param.data = lfs;
param.wt = W; % weighting matrix



% non-linear conjugate gradient method
sus = nlcg(zeros(Nx,Ny,Nz), param);

% OR give simple truncation result as initial guess
% tmp = fftn(lfs)./D;
% T = 0.15; % truncation level
% tmp(abs(D)<T) = 0;
% sus_dc = real(mask.*ifftn(tmp));
% sus = nlcg(sus_dc, param);

sus = real(sus).*mask;
% if want to keep the dipole fitting result
% don't mask it, instead, use the following:
% sus = real(sus);

end


% Rodrigues' rotation formula
function R = rotationmatrix(u,v)
    uvangle=acos(dot(u,v)); %calculate the angle between u and v
    k=cross(u,v); %calculate the cross product between u and v
    k=k./norm(k);
    kx=[0 -k(1,3) k(1,2);k(1,3) 0 -k(1,1);-k(1,2) k(1,1) 0]; %the cross product in matrix form
    R=eye(3,3)+(kx.*sin(uvangle))+((1-cos(uvangle))*kx^2); %Calculate the rotation matrix with Rodrigues' rotation formula
end
