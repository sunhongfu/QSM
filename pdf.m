function lfs = pdf(tfs,mask,par,weight)

[np,nv,ns] = size(tfs);

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

M = 1-mask; % region outside the mask
W = weight.*mask; % weights


b = M.*ifftn(D.*fftn(W.*W.*tfs));
b = b(:);

res = cgs(@Afun,b,1e-6, 200);
lfs = mask.*real(tfs-ifftn(D.*fftn(M.*reshape(res,[np,nv,ns]))));

function y = Afun(x)
    x = reshape(x,[np,nv,ns]);
    y = M.*ifftn(D.*fftn(W.*W.*ifftn(D.*fftn(M.*x))));
    y = y(:);
end

end
