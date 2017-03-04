
function field = forward_field_calc(sus, vox, z_prjs)

%   SUS    : susceptibility distribution
%   VOX    : voxel size, e.g. [1, 1, 1] for isotropic resolution
%   Z_PRJS : normal vector of the imaging plane, e.g. [0, 0, 1] for pure axial


if ~ exist('vox','var') || isempty(vox)
    vox = [1, 1, 1]; % isotropic resolution
end

if ~ exist('z_prjs','var') || isempty(z_prjs)
    z_prjs = [0, 0, 1]; % PURE axial slices
end


imsize = size(sus);

Nx = imsize(1);
Ny = imsize(2);
Nz = imsize(3);

% create K-space filter kernel D
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



field = real(ifftn(D.*fftn(sus)));
