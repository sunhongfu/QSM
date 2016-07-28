% || m * F_{-1} * C * D * F * chi - m * F_{-1} * C * F * B_t|| + alpha * ||chi|| + beta * ||F_{-1} * D * F * chi||

% combined inversion (RESHARP + INV) using CGS
function chi = combined_inversion_cgs(tfs, mask, vox, z_prjs, ker_rad, alpha, beta)


if ~ exist('vox','var') || isempty(vox)
    vox = [1,1,1];
end

if ~ exist('ker_rad','var') || isempty(ker_rad)
    ker_rad = 3;
end

if ~ exist('beta','var') || isempty(beta)
    beta = 1e-3;
end

if ~ exist('cgs_num','var') || isempty(cgs_num)
    cgs_num = 200;
end

imsize = size(tfs);

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
mask_ero = zeros(imsize);
mask_tmp = convn(mask,ker,'same');
mask_ero(mask_tmp > 1-1/sum(h(:))) = 1; % no error points tolerence 


% prepare convolution kernel: delta-ker
dker = -ker;
dker(rx+1,ry+1,rz+1) = 1-ker(rx+1,ry+1,rz+1);
DKER = fftn(dker,imsize); % dker in Fourier domain



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



% b = mask.*ifftn(D.*conj(DKER).*fftn(circshift(mask_ero.*circshift(ifftn(DKER.*fftn(tfs)),-csh),csh)));
b = ifftn(D.*conj(DKER).*fftn(circshift(mask_ero.*circshift(ifftn(DKER.*fftn(tfs)),-csh),csh)));
b = b(:);

% b = H'*(H*tfs(:));
m = cgs(@Afun, b, 1e-6, cgs_num);

chi = real(reshape(m,imsize));

data_fidelity = mask_ero.*circshift(ifftn(DKER.*fftn(tfs-reshape(m,imsize))),-csh);
data_fidelity = norm(data_fidelity(:));

regularization_term = norm(m);

% nested function
function y = Afun(x)
    % y = H'*(H*x) + beta*x;
    x = reshape(x,imsize);

    y = ifftn(D.*conj(DKER).*fftn(circshift(mask_ero.*circshift(ifftn(DKER.*D.*fftn(x)),-csh),csh))) ...
     + alpha.*x +beta.*ifftn(D.*D.*fftn(x));

    y = y(:);
end


end

