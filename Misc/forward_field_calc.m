
function [field, D, d, field_kspace] = forward_field_calc(sus, vox, z_prjs, padding_flag)

%   SUS    : susceptibility distribution
%   VOX    : voxel size, e.g. [1, 1, 1] for isotropic resolution
%   Z_PRJS : normal vector of the imaging plane, e.g. [0, 0, 1] for pure axial


if ~ exist('vox','var') || isempty(vox)
    vox = [1, 1, 1]; % isotropic resolution
end

if ~ exist('z_prjs','var') || isempty(z_prjs)
    z_prjs = [0, 0, 1]; % PURE axial slices
end

if ~ exist('padding_flag','var') || isempty(padding_flag)
    padding_flag = 1; % pad zeros to double size to avoid wrap around
end

imsize = size(sus);

if padding_flag
    sus = padarray(sus, imsize/2);
    imsize = size(sus);
end

Nx = imsize(1);
Ny = imsize(2);
Nz = imsize(3);

% create K-space filter kernel D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) k-space domain
% FOV = vox.*[Nx,Ny,Nz];
% FOVx = FOV(1);
% FOVy = FOV(2);
% FOVz = FOV(3);

% x = -Nx/2:Nx/2-1;
% y = -Ny/2:Ny/2-1;
% z = -Nz/2:Nz/2-1;
% [kx,ky,kz] = ndgrid(x/FOVx,y/FOVy,z/FOVz);
% % D = 1/3 - kz.^2./(kx.^2 + ky.^2 + kz.^2);
% D = 1/3 - (kx.*z_prjs(1)+ky.*z_prjs(2)+kz.*z_prjs(3)).^2./(kx.^2 + ky.^2 + kz.^2);
% D(floor(Nx/2+1),floor(Ny/2+1),floor(Nz/2+1)) = 0;
% D = fftshift(D);

% % (2) image space domain
[X,Y,Z]=ndgrid(-Nx/2:(Nx/2-1),-Ny/2:(Ny/2-1),-Nz/2:(Nz/2-1));

X = X*vox(1);
Y = Y*vox(2);
Z = Z*vox(3);

d = (3*( X*z_prjs(1) + Y*z_prjs(2) + Z*z_prjs(3)).^2 - X.^2-Y.^2-Z.^2)./(4*pi*(X.^2+Y.^2+Z.^2).^2.5);

d(isnan(d)) = 0;
D = fftn(fftshift(d));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

field_kspace = ifftshift(D.*fftn(sus));
field = real(ifftn(D.*fftn(sus)));
D = ifftshift(D);

if padding_flag
    field = field(1+Nx/4:end-Nx/4, 1+Ny/4:end-Ny/4, 1+Nz/4:end-Nz/4);
    D = D(1+Nx/4:end-Nx/4, 1+Ny/4:end-Ny/4, 1+Nz/4:end-Nz/4);
    d = d(1+Nx/4:end-Nx/4, 1+Ny/4:end-Ny/4, 1+Nz/4:end-Nz/4);
    field_kspace = field_kspace(1+Nx/4:end-Nx/4, 1+Ny/4:end-Ny/4, 1+Nz/4:end-Nz/4);
end