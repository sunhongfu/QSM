function [lfs,mask_ero] = pdf(tfs,mask,vox,ker_rad,weight,z_prjs)


if ~ exist('z_prjs','var') || isempty(z_prjs)
    z_prjs = [0, 0, 1]; % PURE axial slices
end


% add zero slices
tfs = padarray(tfs,[0 0 20]);
mask = padarray(mask,[0 0 20]);
weight = padarray(weight,[0 0 20]);



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




M = 1-mask; % region outside the mask
W = weight.*mask; % weights


b = M.*ifftn(D.*fftn(W.*W.*tfs));
b = b(:);

res = cgs(@Afun,b,1e-6, 200);
lfs = mask.*real(tfs-ifftn(D.*fftn(M.*reshape(res,[Nx,Ny,Nz]))));

function y = Afun(x)
    x = reshape(x,[Nx,Ny,Nz]);
    y = M.*ifftn(D.*fftn(W.*W.*ifftn(D.*fftn(M.*x))));
    y = y(:);
end


% erode the edge
rx = round(ker_rad/vox(1));
ry = round(ker_rad/vox(2));
rz = round(ker_rad/vox(3));
[X,Y,Z] = ndgrid(-rx:rx,-ry:ry,-rz:rz);
h = (X.^2/rx^2 + Y.^2/ry^2 + Z.^2/rz^2 < 1);
ker = h/sum(h(:));
csh = [rx,ry,rz];
imsize = size(mask);
mask_tmp = circshift(real(ifftn(fftn(mask).*fftn(ker,imsize))),-csh);
mask_ero = zeros(imsize);
mask_ero(mask_tmp > 1-8/sum(h(:))) = 1; % no error tolerance
lfs = lfs.*mask_ero;


% remove added zero slices
lfs = lfs(:,:,21:end-20);
mask_ero = mask_ero(:,:,21:end-20);

end
