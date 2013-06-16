function [lfs,mask_ero] = pdf(tfs,mask,res,ker_rad,weight)

[np,nv,ns] = size(tfs);


%%%%%%%%%%%%%%%%%%%%% needs scaling with FOV %%%%%%%%%%%%%%%%%%%%%%%
% create K-space filter kernel D
FOV = res.*[np,nv,ns];
FOVx = FOV(1);
FOVy = FOV(2);
FOVz = FOV(3);

x = -np/2:np/2-1;
y = -nv/2:nv/2-1;
z = -ns/2:ns/2-1;
[kx,ky,kz] = ndgrid(x/FOVx,y/FOVy,z/FOVz);
D = 1/3 - kz.^2./(kx.^2 + ky.^2 + kz.^2);

D = fftshift(D);
D(1,1,1) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 1-mask; % region outside the mask
W = weight.*mask; % weights


b = M.*ifftn(D.*fftn(W.*W.*tfs));
b = b(:);

res = cgs(@Afun,b,1e-6, 100);
lfs = mask.*real(tfs-ifftn(D.*fftn(M.*reshape(res,[np,nv,ns]))));

function y = Afun(x)
    x = reshape(x,[np,nv,ns]);
    y = M.*ifftn(D.*fftn(W.*W.*ifftn(D.*fftn(M.*x))));
    y = y(:);
end


% erode the edge
rx = round(ker_rad/res(1));
ry = round(ker_rad/res(2));
rz = round(ker_rad/res(3));
[X,Y,Z] = ndgrid(-rx:rx,-ry:ry,-rz:rz);
h = (X.^2/rx^2 + Y.^2/ry^2 + Z.^2/rz^2 < 1);
ker = h/sum(h(:));
csh = [rx,ry,rz]
imsize = size(mask);
mask_tmp = circshift(real(ifftn(fftn(mask).*fftn(ker,imsize))),-csh);
mask_ero = zeros(imsize);
mask_ero(mask_tmp > 1-6/sum(h(:))) = 1; % 5 voxels error tolerance
lfs = lfs.*mask_ero;



end
